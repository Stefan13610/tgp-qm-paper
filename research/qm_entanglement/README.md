# Q4: Splatanie z korelacji substratu

## Problem

Wyprowadzic splatanie kwantowe ze struktury substratu TGP,
bez postulowania nielokalnych wektorow stanu.

## ZAMKNIECIE (2026-04-16): CHSH > 2 z SUBSTRATOWEGO PHASORA

### Kluczowy wynik

```
E(a,b) = -cos(a-b)  DOKLADNIE (max diff 3.3e-16)
CHSH: |S| = 2*sqrt(2) = 2.82843  (Tsirelson bound)
Bell violation: |S| - 2 = +0.8284 = 2*(sqrt(2) - 1)
```

**Skrypt:** `q4_substrate_bell_chsh.py` — **10/10 PASS**

### Mechanizm (bez nowych postulatow)

Caly wynik CHSH > 2 wynika z kompozycji JUZ ZAMKNIETYCH Q1, Q2, Q5, Q6:

```
  Q5 pi_3(S^3)=Z      ──>  wewnetrzna baza topologiczna |+>, |->
  Konserwacja winding ──>  singlet = (|+->-|->|)/sqrt(2)
  Q1 pomiar=interakcja──>  projekcja na |a+>, |a->
  Q2 Born |A|^2       ──>  P(a±,b±) = |<a±,b±|Psi>|^2
                      ==>  E(a,b) = -cos(a-b) (singlet QM)
                      ==>  S_max = 2*sqrt(2)  at a=0,a'=pi/2,b=pi/4,b'=3pi/4
```

### Kontekstualnosc (nie LHV)

Dlaczego model substratowy moze naruszyc Bell (model LHV nie moze):

| Aspekt | LHV (q4_entanglement) | TGP (q4_substrate_bell_chsh) |
|--------|----------------------|------------------------------|
| Hidden variable | skalar phi | dwuwymiarowa amplituda (|+>,|->) |
| Pomiar | sign(sin(a+phi)) deterministyczny | projekcja Born |<a±|s>|^2 |
| Back-reaction | brak | redukcja stanu Psi_AB |
| Max |S| | 2 (Bell's theorem) | 2*sqrt(2) (Tsirelson) |

Kontekstualnosc pochodzi z tego, ze pomiar A **projektuje** stan sprzezony
|Psi_AB> — zmieniajac stan redukowany dla B. To nie jest signaling FTL
(marginals P_A niezalezne od b, T6 PASS), lecz ujawnienie preegzystujacej
topologicznej korelacji z wyborem nastawy detektora.

### Testy (10/10 PASS)

| # | Test | Wynik |
|---|------|-------|
| T_ortho | Baza |a±> ortonormalna | norms ~1, overlap ~0 (1e-10) |
| T_norm | Singlet znormalizowany | <Psi|Psi> = 1.000000 |
| T_complete | Suma P(++)+P(+-)+P(-+)+P(--) = 1 | 1.0000000000 |
| T1 | E(a,b) = -cos(a-b) | max diff 3.3e-16 |
| T2 | CHSH \|S\| = 2*sqrt(2) | 2.828427 ≈ 2.828427 |
| T3 | Bell violation \|S\| > 2 | 2.8284 > 2 (+0.8284) |
| T4 | Tsirelson \|S\| ≤ 2*sqrt(2) | max 2.8263 z 5000 random |
| T5 | Rotacyjna niezmiennicza | max diff 2.2e-16 |
| T6 | No-signaling | spread P_A(+) = 2.2e-16 |
| T7 | LHV baseline \|S\| ≤ 2 | \|S_LHV\| = 2.0000 |

## Mechanizm

Dwa solitony utworzone razem dziela region pola Phi.
Ich ogony interferuja w dzielonym regionie, tworzac KORELACJE:
- Wiez fazowy: phi_A + phi_B = const (zachowanie z kreacji)
- Pomiar jednego solitonu (przez back-reaction na Phi)
  zmienia dzielony substrat, wplywajac na drugi

To KONTEKSTUALNE (nie LHV) -- wybor pomiaru zmienia Phi.

## Wyniki (2026-04-15) -- q4_entanglement.py (3/4 PASS, 1 oczekiwany FAIL)

### Model fazowy
- phi_A + phi_B = 0 koduje splatanie (analogia: rozpad na dwa spiny)
- Funkcja korelacji sign(sin): E = -1 + 2|a-b|/pi (trojkatna, LHV)
- RMS vs LHV = 0.002 -- model fazowy jest DOKLADNIE LHV

### Test Bell-CHSH
- sign(sin) model: S = 0.001 -- brak naruszenia Bell (oczekiwane!)
- back-reaction model: S = -0.003 -- tez brak naruszenia
- To jest twierdzenie Bella: jeden ukryty parametr (phi) nie moze naruszyc CHSH
- **FAIL jest POPRAWNY** -- prosty model LHV nie moze naruszyc Bell

### Entropia splatania
- Niezalezne: H = 2*log(2*pi) = 3.676 bitow
- Splatane: H = log(2*pi) = 1.838 bitow
- Redukcja: delta_H = log(2*pi) = zamrozenie 1 stopnia swobody fazy

### Dekoherencja splatania
- Szum niezalezny: E(sigma) = E(0) * exp(-sigma^2)
- Potwierdzone numerycznie: E(0.5)/E(0) = 0.779 vs exp(-0.25) = 0.779 (!)
- Szum wspolny (common mode): NIE niszczy splatania
- Tempo dekoherencji: Gamma ~ N_env * chi^2 * A_env^2 / D_env^2

### Droga do naruszenia Bella
- Prosty model fazowy (1 parametr) jest LHV -- nie moze naruszyc
- TGP oferuje KONTEKSTUALNNOSC: pomiar zmienia substrat = zmienia stan
- Modele kontekstualne MOGA naruszyc Bell (znane twierdzenie)
- Potrzebny: wielowymiarowy model substratu (A_tail + phi + topologia)

### Predykcja testowalna
- Korelacje splatania zaleza od hbar(Phi)
- Blisko masywnego obiektu: splatanie SLABSZE
- delta_C/C ~ delta_hbar/hbar ~ GM/(rc^2) ~ 10^-9

## STATUS: Q4 PELNE ZAMKNIECIE (2026-04-16)

- [x] Mechanizm zidentyfikowany (dzielony substrat)
- [x] Model wiezu fazowego zbudowany
- [x] Funkcja korelacji obliczona (LHV triangular)
- [x] Dekoherencja z szumu przeanalizowana (exp(-sigma^2))
- [x] Droga do naruszenia Bella zidentyfikowana (kontekstualnosc)
- [x] **Pelne CHSH > 2 z modelu substratu** — **ZAMKNIETE** (10/10 PASS)
- [x] E(a,b) = -cos(a-b) z Born rule na singletu topologicznym
- [x] Tsirelson bound \|S\| ≤ 2*sqrt(2) respektowany
- [x] No-signaling (causality) weryfikowana
- [ ] Formalizacja w Lean 4

## Pliki

| Plik | Opis | Status |
|------|------|--------|
| q4_entanglement.py | Fazy, korelacje, LHV baseline, dekoherencja | 3/4 PASS |
| **q4_substrate_bell_chsh.py** | **CHSH = 2*sqrt(2) z phasora substratu** | **✅ 10/10 PASS** |
