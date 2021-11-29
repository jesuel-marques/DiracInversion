#!/bin/bash
mkdir MZ Deltas DiracInverses DiracOperators
gcc -o InversaoDiracLivre Fonte/RedeSU3.c Fonte/FuncoesMatematicas.c Fonte/mtwister.c Fonte/SU3.c Fonte/FuncoesInversaoDirac.c InversaoDiracLivre.c Fonte/ranlxd.c Fonte/ranlux_common.c Fonte/MedicoesSU3.c -lm -O3
./InversaoDiracLivre
gcc -o InversaoDiracInteragente Fonte/RedeSU3.c Fonte/FuncoesMatematicas.c Fonte/mtwister.c Fonte/SU3.c Fonte/FuncoesInversaoDirac.c InversaoDiracInteragente.c Fonte/ranlxd.c Fonte/ranlux_common.c Fonte/MedicoesSU3.c -lm -O3
./InversaoDiracInteragente
