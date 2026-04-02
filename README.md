# Topologiczna Detekcja Bifurkacji w Terapii Adaptacyjnej
## Problem Kliniczny
 Terapie nowotworowe nie są idealne. Najczęściej używanym sposobem podawania leku jest tzw. **metoda maksymalnej dawki tolerowanej (MDT)**. Wybijamy jak najwięcej komórek agresywnych, jednak w przypadku guzów heterogennych masowo eliminujemy jedynie komórki wrażliwe na lek. Pozostawia to pustą przestrzeń dla subpopulacji opornej, która szybko proliferuje. Skutkuje to zjawiskiem "uwolnienia konkurencyjnego" (competitive release) i nieuniknionym relapsem.<br>
 <br>
 Częściowym rozwiązaniem jest **terapia adaptacyjna** , zaproponowana w 2009 roku przez R. Gatenby'ego z Moffitt Cancer Center. Zmienia ona cel leczenia: modyfikuje chorobę nowotworową w chorobę przewlekłą. Poprzez regularne cykle ON/OFF dawkowania leku (sterowane stężeniem biomarkera we krwi), terapia utrzymuje konkurencję międzykomórkową i fizycznie ogranicza wzrost guza. Choć mutacje są tu wciąż nieuniknione, pacjenci zyskują cenne lata życia. Celem tego projektu jest użycie topologii algebraicznej do zoptymalizowania detekcji momentu, w którym guz bezpowrotnie mutuje.

## Model Matematyczny: Nieliniowy Układ Hybrydowy
Model zastosowany w tym projekcie opiera się na raku prostaty, ponieważ jest to najlepiej przebadany klinicznie przypadek w kontekście terapii adaptacyjnej. Kluczową rolę odgrywa tutaj marker **PSA** (prostate-specific antigen), którego poziom można regularnie mierzyć we krwi i który stanowi przybliżony wskaźnik wielkości guza.
Użyłam układu czterech sprzężonych równań różniczkowych z publikacji *Zhang et al. (2017)*, reprezentujących dynamikę subpopulacji komórkowych (zależnych od androgenów, niezależnych i produkujących testosteron) oraz dynamikę samego biomarkera (PSA). Zarówno parametry, jak i wartości początkowe zostały zaczerpnięte z oryginalnych materiałów udostępnionych wraz z badaniami klinicznymi [GitHub *Zhang (2017)*](https://github.com/cunninghamjj/Integrating-evolutionary-dynamics-into-treatment-of-mCRPC).

$$
\frac{dT^+}{dt} = r_+ T^+ \left( 1 - \frac{T^+ + \alpha_{12} T^P + \alpha_{13} T^-}{K_+(t)} \right)
$$

$$
\frac{dT^P}{dt} = r_P T^P \left( 1 - \frac{\alpha_{21} T^+ + T^P + \alpha_{23} T^-}{K_P} \right)
$$

$$
\frac{dT^-}{dt} = r_- T^- \left( 1 - \frac{\alpha_{31} T^+ + \alpha_{32} T^P + T^-}{K_-} \right)
$$

$$
\frac{d(PSA)}{dt} = (T^+ + T^P + T^-) - \sigma \cdot PSA
$$

Gdzie poszczególne zmienne oznaczają:<br>
$T^+$, $T^P$, $T^-$ – populacje komórek: zależnych od androgenów, produkujących testosteron oraz niezależnych (lekoopornych).<br>
$r_i$ – współczynniki tempa proliferacji dla poszczególnych typów komórek.<br>
$K_i$ – pojemność środowiska. Co istotne, $K_+(t)$ jest zmienną zależną od czasu - jej wartość ulega skokowej zmianie w zależności od tego, czy terapia jest włączona, czy wyłączona, co tworzy z układu system hybrydowy.<br>
$\alpha_{ij}$ – parametry macierzy Lotki-Volterry określające presję konkurencyjną, jaką populacja $j$ wywiera na populację $i$<br>
$\sigma$ – stała szybkości rozpadu biomarkera PSA.<br>

<p> To nieliniowy, przełączany układ hybrydowy. Sterowanie dawką leku wprowadza skokowe zmiany pojemności środowiska $K$ w równaniach Lotki-Volterry. Z numerycznego punktu widzenia to układ silnie sztywny (stiff), charakteryzujący się wieloskalowością czasu (wolny wzrost komórek kontra szybki rozpad PSA), co wymusiło użycie stabilnych, niejawnych metod całkowania wstecznego (BDF) </p>

## Metodologia: Topologiczna Analiza Danych (TDA)
### Twierdzenie Takensa (Zanurzenie z opóźnieniem)
Twierdzenie Takensa, udowodnione w latach 80., przedstawia zestaw założeń, które umożliwiają nam rekonstrukcję złożonego układu dynamicznego za pomocą wyłącznie jednej jego obserwowalnej współrzędnej (zależnej od czasu). U mnie konkretnie umożliwia ono odtworzenie wielowymiarowej dynamiki całego modelu (czyli ukrytych wewnątrz guza populacji komórek) patrząc wyłącznie na szereg czasowy zmian biomarkera z krwi.
<br>
Zanurzając te wartości wraz z opóźnieniami w odpowiednio wymiarowej przestrzeni, dostajemy matematyczną gwarancję, że nowy układ będzie miał te same niezmienniki topologiczne, co układ oryginalny. Osiągamy to poprzez konstrukcję wektora zanurzenia (embedding vector): 

$$
\mathbf{x}(t) = [PSA(t), PSA(t+\tau), PSA(t+2\tau), \dots, PSA(t+(m-1)\tau)]
$$

gdzie $\tau$ to precyzyjnie dobrane opóźnienie czasowe (oparte na średnim okresie zdrowego cyklu), a $m$ to wymiar przestrzeni fazowej.
### Homologia Persystentna
Homologia persystentna jest kolejnym narzędziem umożliwiającym analizę naszego zrekonstruowanego układu. Zamiast skupiać się na konkretnych, dokładnych wartościach liczbowych, analizuję portret fazowy pod kątem występowania i struktury cykli (pętli).<br>
Dla otrzymanej chmury punktów w przestrzeni fazowej konstruuję tzw. kompleksy symplicjalne, wykorzystując filtrację Vietorisa-Ripsa. Śledzimy ewolucję pierwszej liczby Bettiego ($\beta_1$), która w topologii zlicza jednowymiarowe "dziury". W kontekście biologicznym, taka stabilna, topologiczna dziura odpowiada prawidłowemu cyklowi leczenia pacjenta.<br>

Wysoka wartość persystencji dla $\beta_1$ daje nam kategoryczną odpowiedź, że pacjent poprawnie reaguje na terapię adaptacyjną, a atraktor posiada zamkniętą geometrię. Spadek tej wartości to matematyczny dowód na przerwanie cyklu.<br>
## Algorytm Detekcji: "Sliding Window"
Liczenie homologii z portretu fazowego całego układu nie dostarcza informacji o zmianach w czasie. Dlatego kluczowym narzędziem informatycznym w projekcie jest algorytm okna przesuwnego (sliding window).<br>
Tnę szereg czasowy PSA na odcinki i za pomocą okna o określonej szerokości przechodzę przez sygnał, zapisując maksymalną żywotność (persystencję) jednowymiarowych cykli. Tworząc z tych wartości wykres w czasie, otrzymuję dowód na wystąpienie bifurkacji układu. Detektor topologiczny jednoznacznie wskazuje moment załamania się cyklu granicznego – punkt, w którym guz przestaje reagować na terapię i następuje ucieczka nowotworu. Eliminuje to potrzebę subiektywnej oceny wizualnej wykresów.<br>

## Wyniki Działania Algorytmu

## Perspektywy i Ograniczenia
W porównaniu do standardowych algorytmów analizy sygnałów (jak np. autokorelacja), topologia ma gigantyczny atut: nie zależy od dokładnych wartości i jest odporna na deformacje układu. Omija problem żmudnego doboru parametrów, który do dziś stanowi wyzwanie w literaturze przedmiotu.
<br>
Głównym ograniczeniem TDA jest zapotrzebowanie na gęste dane. Obecnie pomiary krwi wykonuje się co kilka tygodni lub miesięcy, co nie pozwala na pełną rekonstrukcję przestrzeni fazowej. Przyszłym rozwiązaniem dla wdrożenia tego detektora mogą być zaawansowane metody numeryczne (inteligentna interpolacja) lub rozwój biofizycznych urządzeń typu wearables, zdolnych do ciągłego, nieinwazyjnego monitorowania biomarkerów.
## Literatura
1. Zhang, J. Z., Cunningham, J. J., Brown, J. S., & Gatenby, R. A. (2017). *Integrating evolutionary dynamics into treatment of
metastatic castrate-resistant prostate cancer*.
2. Gatenby, R. A. (2009). *Adaptive Therapy*.
3. Salvioli, M., Vandelaer, L., Baena, E., Schneider, K., Cavill, R., & Staňková, K.  
*The effect of tumor composition on the success of adaptive therapy: The case of metastatic castrate-resistant prostate cancer*.
