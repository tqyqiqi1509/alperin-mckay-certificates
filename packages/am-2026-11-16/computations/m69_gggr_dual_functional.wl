em[i_,j_]:=Normal[SparseArray[{{i,j}->1},{8,8}]];
xm[i_,j_]:=em[i,j]-em[j+4,i+4];
xp[i_,j_]:=em[i,j+4]-em[j,i+4];
xx=xm[1,2]+xm[2,3]+xp[2,3];
xd=-Transpose[xx];
hh=DiagonalMatrix[{4,2,0,0,-4,-2,0,0}];
gg=IdentityMatrix[8][[All,{1,2,3,8,5,6,7,4}]];
basis2={xm[1,2],xm[2,3],xp[2,3],xm[2,4],xp[2,4]};
coef=(Tr[xd.#]& /@ basis2);
out=<|"xGradeTwo"->(hh.xx-xx.hh==2xx),"dualGradeMinusTwo"->(hh.xd-xd.hh==-2xd),"graphFixesX"->(gg.xx.gg==xx),"graphFixesDual"->(gg.xd.gg==xd),"graphFixesGrading"->(gg.hh.gg==hh),"traceFunctionalCoefficients"->coef,"coefficientGCD"->Apply[GCD,coef],"positivePairingCoefficients"->(Tr[xx.#]& /@ basis2),"selfPairing"->Tr[xd.xx],"functionalStillNonzeroModuloThree"->(Mod[coef,3]!={0,0,0,0,0}),"rootOrbitPermutation"->(First[FirstPosition[basis2,gg.#.gg]]& /@ basis2)|>;
ExportString[out,"RawJSON"]
