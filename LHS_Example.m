%% ================================================================
%% Basit meteoroloji örneği: 4 değişkenli LHS 
%% Değişkenler (girdiler): 
%%   - T(°C)  : Sıcaklık
%%   - p(hPa) : Basınç
%%   - RH(%)  : Bağıl nem (Relative Humidity)
%%   - U(m/s) : Rüzgâr hızı
%%
%% Çıktılar (hesaplananlar):
%%   - rho (kg/m^3): Hava yoğunluğu (nem etkisi dahil)
%%   - a (m/s)     : Ses hızı (kuru hava yaklaşımıyla)
%%   - q (Pa)      : Dinamik basınç = 0.5 * rho * U^2
%%
%% Notlar:
%%  - ".*" eleman-eleman çarpma demektir (vektör/matrislerde kullanılır).
%%  - "./" eleman-eleman bölme.
%%  - ".^" eleman-eleman üstel alma (kuvvet).
%%  - Tek tırnak (') burada transpoz (satır<->sütun) için kullanılıyor; 
%%    sayılar reel olduğu için kompleks eşlenik konusu önemli değil.
%% ================================================================

rng(42);               % Rastgele sayı üreticisini sabitliyoruz.
                       % Böylece her çalıştırışta aynı "rastgele" sonuçları alırız (tekrar üretilebilirlik).

n = 500;               % n = kaç tane örnek (satır) üreteceğiz? (örnek sayısı)
d = 4;                 % d = kaç boyut/kaç değişken var? (T, p, RH, U) -> 4

% LHS (Latin Hypercube Sampling) ile [0,1]^d içinde n adet örnek üret.
% latin_hypercube fonksiyonunu aşağıda tanımlıyoruz; toolbok gerektirmez.

Ulhs = latin_hypercube(n, d);   

% -----------------------------
% 1) GİRİŞ DEĞİŞKENLERİNİN ARALIKLARI
% -----------------------------
% Burada her değişkenin MIN-MAX sınırını tanımlıyoruz (basit, gerçekçi değerler).

T_min  = 0;    T_max  = 40;     % Sıcaklık aralığı [0, 40] °C
p_min  = 980;  p_max  = 1030;   % Basınç aralığı [980, 1030] hPa (deniz seviyesine yakın)
RH_min = 10;   RH_max = 90;     % Bağıl nem [%] aralığı [10, 90]
U_min  = 0;    U_max  = 15;     % Rüzgâr hızı [0, 15] m/s (0: sakin, 15: kuvvetli esinti)

% -----------------------------
% 2) [0,1] ARALIĞINI GERÇEK DEĞERLERE ÖLÇEKLEME
% -----------------------------
% Ulhs(:,k) şu an 0-1 arasında değerler. Bunları kendi fiziksel aralıklarına taşıyoruz.
% Formül: fiziğe_çıktı = Ulhs * (MAX - MIN) + MIN

T_C   = Ulhs(:,1) * (T_max  - T_min ) + T_min ;   % °C cinsinden sıcaklık vektörü (n x 1)
p_hPa = Ulhs(:,2) * (p_max  - p_min ) + p_min ;   % hPa cinsinden basınç vektörü
RH    = Ulhs(:,3) * (RH_max - RH_min) + RH_min;   % % cinsinden bağıl nem vektörü
U_ms  = Ulhs(:,4) * (U_max  - U_min ) + U_min ;   % m/s cinsinden rüzgâr hızı vektörü

% -----------------------------
% 3) BİRİM DÖNÜŞÜMLERİ
% -----------------------------

T_K  = T_C + 273.15;      % °C -> K çeviriyoruz (Kelvin = Celsius + 273.15)
p_Pa = p_hPa * 100;       % hPa -> Pa çeviriyoruz (1 hPa = 100 Pa)

% -----------------------------
% 4) SABİTLER
% -----------------------------

Rd    = 287.05;   % Kuru hava için gaz sabiti (J/(kg*K))
Rv    = 461.495;  % Su buharı için gaz sabiti (J/(kg*K))
gamma = 1.4;      % Kuru hava için Cp/Cv (ısı oranı) ~1.4 kabul edilir.

% -----------------------------
% 5) NEM HESABI (Magnus-Tetens yaklaştırımı)
% -----------------------------
% Amaç: Bağıl nem (RH) ve sıcaklıktan (T_C) su buharı kısmi basıncını bulmak.
% Adım 1: Doygunluk buhar basıncı (e_s) [hPa] -> sıcaklığa bağlıdır.
% Formül: e_s = 6.112 * exp( (17.62*T) / (243.12 + T) )  -> T: °C
% Geçerlilik: Yaklaşık -45°C ile +60°C aralığında iyi çalışır.

e_s_hPa = 6.112 * exp((17.62 .* T_C) ./ (243.12 + T_C)); 

% Adım 2: Gerçek buhar basıncı e [hPa] = (RH/100) * e_s
% Yani RH=100% ise hava buhar açısından doygun; RH düşükse e daha küçük.
e_hPa = (RH / 100) .* e_s_hPa;  

% Pa cinsine geçiyoruz çünkü ana basıncı p_Pa (Pa) cinsinden kullanacağız.

e_Pa = e_hPa * 100;        

% -----------------------------
% 6) HAVA YOĞUNLUĞU (rho) HESABI
% -----------------------------
% İki gazın (kuru hava + su buharı) karışımı gibi düşün.
% Toplam basınç: p = p_d + e  (p_d: kuru havanın basıncı, e: su buharının basıncı)
% Yoğunluk: rho = (p_d)/(Rd*T) + (e)/(Rv*T)
% Burada T Kelvin cinsinden, p_d ve e Pascal (Pa) cinsinden.

rho = ((p_Pa - e_Pa) ./ (Rd .* T_K)) + (e_Pa ./ (Rv .* T_K));  % Sonuç: kg/m^3

% -----------------------------
% 7) SES HIZI (a) HESABI
% -----------------------------
% Basit/kabaca: a = sqrt(gamma * Rd * T)  -> Kuru hava kabulü.
% Not: Nem etkisini hesaba katan daha detaylı formüller var; gerekirse eklenir.
a = sqrt(gamma * Rd .* T_K);   % Sonuç: m/s

% -----------------------------
% 8) DİNAMİK BASINÇ (q) HESABI
% -----------------------------
% Akışkanlar mekaniği: q = 0.5 * rho * U^2
% U burada hız (m/s), rho yoğunluk (kg/m^3), q ise Pascal (Pa).

q = 0.5 .* rho .* (U_ms .^ 2); % ".^2" = eleman eleman karesini al.

% -----------------------------
% 9) KISA ÖZET BASTIRMA
% -----------------------------
% median(...) vektörün ortanca değerini verir.

fprintf('rho medyan = %.3f kg/m^3, a medyan = %.1f m/s, q medyan = %.1f Pa\n', ...
    median(rho), median(a), median(q));

% -----------------------------
% 10) HIZLI GÖRSELLEŞTİRMELER
% -----------------------------
% Histogram: rho dağılımını görmek için 30 aralıklı histogram çiziyoruz.

figure; 
histogram(rho, 30); 
xlabel('\rho (kg/m^3)'); ylabel('Sıklık'); title('Hava yoğunluğu (LHS)');

% Ses hızı dağılımı

figure; 
histogram(a, 30); 
xlabel('a (m/s)'); ylabel('Sıklık'); title('Ses hızı (LHS)');

% Sıcaklık vs Yoğunluk saçılım grafiği (sıcaklık arttıkça yoğunluğun azalma eğilimini görürsün)

figure; 
scatter(T_C, rho, 12, 'filled'); grid on; 
xlabel('T (°C)'); ylabel('\rho (kg/m^3)'); title('T vs. \rho');



%% ================================================================
%% LHS ÜRETİCİ FONKSİYON: latin_hypercube (toolbox gerektirmez)
%% Amaç: [0,1]^d içinde n adet "Latin Hiperküp" örneği üretmek.
%% LHS mantığı:
%%   1) Her ekseni n adet "dilime" ayır (eşit aralıklı).
%%   2) Her dilimden tam 1 nokta seç (o dilimin içinde rastgele bir yer).
%%   3) Bu n noktayı satırlara rastgele sırayla yerleştir (permütasyon),
%%      böylece sütunlar arası yapay hizalanma/korelasyon kırılır.
%% Sonuç: Her sütunda (her değişkende) tüm aralık eşit şekilde kapsanmış olur
%% (her dilim tam 1 kere temsil edilir).
%% ================================================================

function X = latin_hypercube(n, d)
    % Girdi:
    %   n : üreteceğimiz örnek sayısı (satır sayısı)
    %   d : değişken/boyut sayısı (sütun sayısı)
    % Çıktı:
    %   X : [n x d] boyutunda örnek matrisi. 
    %       X(i,j) -> i'inci örneğin j'inci değişken değeri (0..1 arasında)

    % 0'larla doldurulmuş bir [n x d] matris oluşturuyoruz (performans için ön ayırma).

    X = zeros(n, d);  
    
    % edges = [0, 1/n, 2/n, ..., 1] gibi; yani [0,1] aralığını n eşit dilime bölen sınırlar.
    % edges uzunluğu (n+1) olur. Dilim sayısı n olduğu için, sınır sayısı n+1'dir.

    edges = linspace(0, 1, n+1); 

    % Her sütunu (her değişkeni) tek tek kuracağız:

    for j = 1:d
        % u = [n x 1] vektör; 0 ile 1 arasında rastgele sayılar (her dilim için 1 tane).
        % Bu değerler, dilim içindeki "konumu" belirler (tam sol, ortalar, sağ yakın...).

        u = rand(n, 1);  
        
        % AŞAĞIDAKİ SATIR, "HER DİLİMDEN BİR NOKTA" ÜRETİR:
        % edges(1:end-1)  -> [L1, L2, ..., Ln]  = her dilimin sol sınırı
        % diff(edges)     -> [W1, W2, ..., Wn]  = her dilimin genişliği (Wk = edges(k+1)-edges(k))
        % u .* diff(edges)' -> [u1*W1; u2*W2; ...; un*Wn] = her dilim için içeri doğru ofset
        % toplarsak: pts(k) = Lk + u_k * Wk   -> k'inci dilimin içinde rasgele bir nokta
        %
        % Not: Sonuna koyduğumuz (') transpoz demek; satır/sütun boyutlarını uyumlu
        %      hale getirmek için kullanıyoruz (kolon vektör olsun diye).

        pts = edges(1:end-1)' + u .* diff(edges)';  

        % ŞİMDİ PERMÜTASYON:
        % randperm(n) -> 1..n sayılarını rastgele karıştırır (tekrarsız).
        % pts(randperm(n)) -> pts içindeki n noktayı rastgele bir sırada döndürür.
        % Neden? Eğer hep aynı sıra ile koyarsak sütunlar arası hizalanma olur (kötü).
        % Permütasyonla her sütunda dilimler farklı sırada satırlara dağılır -> Latin özelliği korunur,
        % boyutlar arası yapay korelasyon kırılır.

        X(:, j) = pts(randperm(n)); 
    end
end
