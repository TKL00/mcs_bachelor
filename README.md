# Bachelorprojekt - Største Fælles Udviddelse af "Anker"
af Kasper Halkjær Beider og Tobias Klink Lehn, datalogistuderende ved Institut for Matematik og Datalogi på Syddansk Universitet, Odense.

Biblioteket implementerer to forskellige algoritmer til at finde største fælles delgraf (udviddelse). Den ene er en naiv implementering af James McGregor's backtracking-algoritme (/src/mcgregor.py), den anden er en udviddelse af en algoritmer, som bruger produktgrafer/klikefinding ud fra et anker. Sidstnævnte (src/cliques.py) er formuleret af Akbar Davoodi (SDU).

Algoritmen er tiltænkt at finde en udviddelse af en på forhånd kendt kontekst af kemiske molekyler (reaktionskerner). Eksempler på disse molekyler kan findes i _labelled\_graphs_.
