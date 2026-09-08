# Ueberarbeiteter Mailentwurf an Heinrich Krobath

**Warum ueberarbeitet:** Zwei der vier urspruenglichen Fragen sind aus den gelieferten Dateien
beantwortbar, und die Praemisse von Frage 1 ("Frame 300 = 30 ns") ist widerlegt. Wer die alte Mail
abschickt, fragt nach Dingen, die er selbst messen kann, und uebersieht die eigentlich kritische
Frage: **ob es ueberhaupt einen Produktionslauf gibt.**

Belege zu jedem Punkt: `Scripts/Analysis/MD_FRAME_PROVENANCE_ANSWERS_2026-09-03.md`.

---

## Was sich gegenueber dem ersten Entwurf aendert

| urspruenglich | jetzt |
|---|---|
| "Wie lang war der Production-Run?" als offene Frage | Aus den Dateien gemessen: `bar3` umfasst rund 0.33 ns bei ~0.66 ps pro Frame. Die Frage lautet jetzt, **ob `bar3` ueberhaupt der Produktionslauf ist.** |
| "Frame 300 = 30 ns" als Formulierungsziel | Faellt weg. Frame 300 liegt bei rund 0.2 ns. |
| "300, 400, 499 sehen nach systematischer Auswahl aus" | Bestaetigt: 499 **ist** der letzte Frame der Datei. Die Frage nach RNG/Seed eruebrigt sich. |
| "ist Najjar et al. 2025 die richtige Quelle?" | Nein. Das dort beschriebene System ist ein anderes (siehe unten). Die Frage lautet jetzt, welches System dieses ist. |

---

## Entwurf

> **Betreff: Rueckfragen zu den MD-Frames (Fr0/300/400/499) — mit Zwischenstand meinerseits**
>
> Hallo Heinrich,
>
> vielen Dank fuer die Dateien, das hilft mir sehr. Ich habe die DCDs inzwischen selbst ausgelesen
> und komme bei zwei Punkten zu einem Ergebnis, das ich gern mit dir abgleichen wuerde, bevor ich
> es in die Methodik schreibe.
>
> **Was ich messen konnte.** Die DCD-Header sind offenbar mit VMD neu geschrieben worden und
> enthalten nur Platzhalter (`nsavc = 1`, `delta = 1.0`), also keine Zeitinformation. Ich habe das
> Frame-Intervall deshalb ueber die mittlere quadratische Verschiebung der Wassermolekuele
> abgeschaetzt. Danach liegen die Frames in `Orai1_WT_bar3.dcd` rund 0.66 ps auseinander, die Datei
> umfasst also insgesamt etwa 0.33 ns. `Orai1_WT_fix1.dcd` sieht nach einem Aufheizlauf mit
> fixiertem Protein bei konstantem Volumen aus, ungefaehr 38 ps lang.
>
> **Meine eigentliche Frage ist deshalb eine andere geworden:**
>
> 1. **Ist `bar3` der Produktionslauf, oder eine Equilibrierungsstufe?** Der Name und die Dichte
>    der gespeicherten Frames sprechen fuer Letzteres, und zwischen `fix1` und `bar3` fehlt
>    offensichtlich mindestens ein Zwischenschritt: das Protein bewegt sich dazwischen um 5.8 A
>    und die Box schrumpft von 2273 auf 1750 nm3. Gibt es einen laengeren Produktionslauf, aus dem
>    die Snapshots eigentlich stammen sollten? Falls ja, haette ich den gern, weil meine vier
>    Rezeptorstrukturen sonst nur rund 0.33 ns abdecken.
>
> 2. **Zeitschritt und `dcdfreq`.** Wenn du die NAMD-Inputs (`step6.*`, `step7*.inp` oder eure
>    eigenen `fix*.conf` / `bar*.conf`) noch hast, waere mir mit denen am meisten geholfen. Dann
>    kann ich die Zeitangaben belegen statt sie zu schaetzen.
>
> 3. **Homologiemodell — bitte nur kurz bestaetigen.** Ich habe das Modell aus dem ModelArchive
>    (`ma-akdjp`, zu Frischauf et al., Sci. Signal. 8, ra131, 2015) heruntergeladen und auf Fr0
>    gelegt: sechs Ketten, Reste 66-288, keine einzige abweichende Aminosaeure, 1338 C-alpha mit
>    0.12 A RMSD, und die leichte Abweichung von der sechszaehligen Symmetrie am M4-Ende ist in
>    beiden identisch. Ich wuerde also Frischauf 2015 als Primaerquelle zitieren und Hou et al. 2012
>    nur als Template. Passt das so, oder wurde das Modell fuer dieses System noch einmal
>    ueberarbeitet?
>
> 4. **Simulations-Setup.** Najjar et al. 2025 passt als Referenz fuer dieses System nicht: dort
>    ist die Membran DDPC/DLPE/DLPG/DLPS/DMPI/PSM ohne Cholesterin, mit CHARMM-GUI V2.0 von 2020
>    gebaut, 304.917 Atome, Reste 62-297. Das mir gelieferte System hat dagegen eine
>    cholesterinreiche Membran (213 Cholesterin, 71 Sphingomyelin, 71 DLPC, 71 DLPE, also 50 mol%
>    Cholesterin), wurde am 6. April 2022 mit CHARMM-GUI v3.7 gebaut, hat 182.692 Atome, die Reste
>    66-288 und enthaelt nur KCl, kein Kalzium. Gibt es zu **diesem** Aufbau eine Publikation oder
>    zumindest eine Notiz, aus der Kraftfeld, Membranzusammensetzung und Parameter hervorgehen?
>
> Zwei Dinge, die ich nicht mehr fragen muss und die ich so schreiben wuerde, wenn du nicht
> widersprichst: Fr499 ist der letzte Frame der Datei, und Fr300/Fr400/Fr499 liegen in gleichen
> Abstaenden davor, ich beschreibe die Auswahl also als systematisch und nicht als Zufallsziehung.
> Fr0 stammt aus `fix1`, ist also die minimierte Ausgangsstruktur vor Beginn der Dynamik und kein
> Frame des spaeteren Laufs.
>
> Wenn es einfacher ist, klaeren wir das gern kurz telefonisch.
>
> Beste Gruesse
> Dominik

---

## Hinweis zur Umlaut-Schreibweise

Der Entwurf oben ist bewusst ohne Umlaute gesetzt. Vor dem Versenden ae/oe/ue/ss durch ä/ö/ü/ß
ersetzen.
