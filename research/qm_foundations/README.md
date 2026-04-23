# Q0: Emergentna Mechanika Kwantowa z TGP -- Architektura

## Teza centralna

QM nie jest fundamentalna. Wynika z **samozwrotnosci pola Phi**:
czastka (soliton) tworzy osrodek (Phi), w ktorym jest mierzona,
a obserwator (inny soliton) tworzy wlasny osrodek, ktory interferuje z mierzonym.

## Lancuch derywacji

```
WARSTWA 0: Ontologia
  Phi konstytuuje przestrzen (Aksjomat A1)
  Czastka = soliton Phi (sekcja 5-6)
  Pomiar = interakcja soliton-soliton

WARSTWA 1: Niepewnosc pomiarowa (Q1)
  Soliton tworzy metryke g_ij = g*delta_ij
  Pomiar wymaga drugiego solitonu (detektor)
  Nakladanie pol Phi zmienia to, co mierzymy
  => Delta_x * Delta_p >= hbar(Phi)

WARSTWA 2: Regula Borna (Q2)
  Ogon solitonu: g(r) ~ 1 + A*sin(r+delta)/r
  Interferencja ogonow daje sygnal ~ |A|^2
  Powtorzony pomiar => rozklad ~ |A_tail|^2
  => P(x) = |psi(x)|^2

WARSTWA 3: Superpozycja (Q3)
  Rownanie liniowe w przyblizeniu slabego pola
  g = 1 + eps*f => f'' + (2/r)f' + f = 0 (liniowe!)
  Superpozycja = liniowa kombinacja rozwiazan f

WARSTWA 4: Splatanie (Q4)
  Dwa solitony dziela krawedzie substratu Gamma
  Korelacje nielokalne z topologii grafu
  Czy naruszone nierownosci Bella?

WARSTWA 5: Spin, statystyka, dekoherencja (Q5-Q7)
  Spin 1/2 z topologii solitonu (pi_3(S^3))
  Fermi-Dirac vs Bose-Einstein z permutacji
  Dekoherencja z tlumienia hbar(Phi) w gestej przestrzeni
```

## Kluczowa roznica TGP vs standardowa QM

| Aspekt | Standardowa QM | TGP |
|--------|---------------|-----|
| Przestrzen | Tlo (dana) | Dynamiczna (Phi) |
| Czastka | Punkt/fala na tle | Soliton tworzacy tlo |
| Pomiar | Postulat (kolaps) | Interakcja soliton-soliton |
| hbar | Stala uniwersalna | hbar = pi*chi*A_0/sqrt(alpha) (z symetrii skalowania ODE) |
| Born | Postulat |psi|^2 | Emergentna z interferencji ogonow |
| Superpozycja | Aksjomatyczna | Liniownosc rownania w slabym polu |

## Podfoldery badawcze

| ID | Folder | Problem | Priorytet |
|----|--------|---------|-----------|
| Q1 | `qm_measurement/` | Niepewnosc pomiarowa z samozwrotnosci | NATYCHMIAST |
| Q2 | `qm_born_rule/` | Regula Borna z |A_tail|^2 | WYSOKI |
| Q3 | `qm_superposition/` | Liniownosc z perturbacji | SREDNI |
| Q4 | `qm_entanglement/` | Splatanie z substratu | SREDNI |
| Q5 | `qm_spin/` | Spin 1/2 z topologii | NISKI |
| Q6 | `qm_statistics/` | Statystyka czastek | NISKI |
| Q7 | `qm_decoherence/` | Dekoherencja z hbar(Phi) | NISKI |

## Status

- [x] Q1: Niepewnosc pomiarowa — **ZAMKNIETE** (22/25 PASS, 4 skrypty)
- [x] Q2: Regula Borna — **ZAMKNIETE** (wynika z Q1: p=2.028, CV=2.3%)
- [x] Q3: Superpozycja — **ZAMKNIETE** (7/7 PASS, linearyzacja + NL korekty)
- [x] Q4: Splatanie — **ZAMKNIETE** (2026-04-16: q4_substrate_bell_chsh 10/10, CHSH=2*sqrt(2))
- [x] Q5: Spin — **ZAMKNIETE** (7/7 PASS, pi_3(S^3)=Z, B=1, FR -> spin 1/2)
- [x] Q6: Statystyka — **ZAMKNIETE** (8/8 PASS, FD/BE z topologii, anyony w 2D)
- [x] Q7: Dekoherencja — **ZAMKNIETE** (8/8 PASS, 3 drogi, darwinizm kwantowy)

## Wyniki analityczne (q0_analytical.py, 5/5 PASS)

Kluczowe wyniki wyprowadzone ALGEBRAICZNIE z ODE substratu:

1. **Born p=2 ANALITYCZNIE**: Z liniowosci perturbacji chi NIE zalezy od A_part
   => Signal ~ chi^2 * A_part^2 => p=2 dokladnie. Numerycznie: slope=1.014, R^2=0.99998.
   Odchylenie p=2.028 z wyzszych rzedow perturbacji.

2. **chi ~ 0.86 (prawie uniwersalne)**: CV=7% w zakresie g0=[0.3, 0.95].
   chi jest wlasciwoscia RDZENIA solitonu, nie ogona.

3. **hbar = pi*chi*A_tail**: WYPROWADZONE z parametrow ODE, nie postulowane.
   Dx*Dp = pi*chi*A/D -- k kasuje sie w iloczynie niepewnosci.
   SYMETRIA SKALOWANIA ODE: g(r;alpha) = g_ref(sqrt(alpha)*r) [DOKLADNA]
   => A(alpha) = A_0/sqrt(alpha) [wykl. -0.49999, R^2=0.9999999955]
   => hbar(alpha) = hbar_0/sqrt(alpha) [DOKLADNE]
   => hbar(Phi) = hbar_0*sqrt(Phi_0/Phi) jest WYPROWADZONE (nie postulowane!)
   Mechanizm: gestszy substrat => wieksze alpha => mniejsze A => mniejsze hbar.

4. **C_NL = 0.535 (I rzad perturbacji)**: Gorne ograniczenie na korekty NL.
   Pelne ODE daje C_full ~ 0.37 (wyzsze rzedy redukuja o 31%).

5. **Druga harmoniczna A_2/A_1 = 1.1%**: Testowalna sygnatura nieliniowosci.

## Pliki

| Plik | Opis | Status |
|------|------|--------|
| q0_analytical.py | Analityczne derywacje: Born, chi, hbar, C_NL | 5/5 PASS |
| q0_hbar_scaling.py | Symetria skalowania ODE: A~1/sqrt(alpha), hbar(Phi) | 7/7 PASS |
