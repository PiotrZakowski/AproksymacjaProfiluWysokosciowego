/*************************************************************
 *	Copyright (c) 2017 Piotr Pawe³ ¯akowski.
 *	All rights reserved.
 *************************************************************/

#include "InterpolacjaLagrange.h"
#include "InterpolacjaFunkcjami3Stopnia.h"

using namespace std;

/*
 * Struktura zawierajaca nazwy dla plikow majacych zawierac wyniki interpolacji.
 */
struct NazwyPlikowZapisu
{
	char *nazwa64_interpolacjaLagrange;
	char *nazwa64_interpolacjaFS3S;
	char *nazwa16_interpolacjaLagrange;
	char *nazwa16_interpolacjaFS3S;
};

/*
 * Wczytuje plik z danymi i zapamietuje dane okreslonych punktow.
 */
void wyciagnijDaneZPliku_daneOWezlach(DaneFunkcji **daneWezlow, char *nazwaPliku, int liczbaWierszy, int coKtoryWezelBrac)
{
	FILE *plik = fopen(nazwaPliku, "r");

	//usuniecie naglowkow tabel
	char naglowek[20];
	fscanf(plik, "%s (m),%s (m)\n", naglowek, naglowek);

	*daneWezlow = (DaneFunkcji *)malloc((((liczbaWierszy) / coKtoryWezelBrac) + 1) * sizeof(DaneFunkcji));
	double pobranyX = NULL;
	double pobranaFx = NULL;
	int iteratorPliku = -1;
	int iteratorDanych = -1;
	while (!feof(plik))
	{
		iteratorPliku++;
		fscanf(plik, "%lf,%lf\n", &pobranyX, &pobranaFx);
		if (iteratorPliku%coKtoryWezelBrac == 0)
		{
			iteratorDanych++;
			(*daneWezlow)[iteratorDanych].x = pobranyX;
			(*daneWezlow)[iteratorDanych].fx = pobranaFx;
		}

	}
	(*daneWezlow)[((liczbaWierszy) / coKtoryWezelBrac)].x = pobranyX;
	(*daneWezlow)[((liczbaWierszy) / coKtoryWezelBrac)].fx = pobranaFx;

	fclose(plik);

	return;
}

/*
 * Funkcja tworzy profile wysokosciowe za pomoca interpolacji
 */
void wyznaczDaneProfile(char *nazwaPlikuDanych, NazwyPlikowZapisu nazwyPlikuZapisu, double xPoczatkowe, double xKoncowe, double xSkok)
{
	int liczbaRekordow = 512;
	DaneFunkcji *daneWezlow = NULL;
	//wszytskie punkty
	wyciagnijDaneZPliku_daneOWezlach(&daneWezlow, nazwaPlikuDanych, liczbaRekordow, 1);

	//co 64 punkt
	DaneFunkcji *daneWezlow_64 = NULL;
	wyciagnijDaneZPliku_daneOWezlach(&daneWezlow_64, nazwaPlikuDanych, liczbaRekordow, 64);

	iterujPoFunkcji_interpolacjaLagrange(daneWezlow_64, liczbaRekordow / 64, nazwyPlikuZapisu.nazwa64_interpolacjaLagrange, xPoczatkowe, xKoncowe, xSkok);
	szacowanieBleduAproksymacji_interpolacjaLagrange(daneWezlow_64, liczbaRekordow / 64, nazwyPlikuZapisu.nazwa64_interpolacjaLagrange, daneWezlow, liczbaRekordow);

	iterujPoFunkcji_interpolacjaFS3S(daneWezlow_64, liczbaRekordow / 64, nazwyPlikuZapisu.nazwa64_interpolacjaFS3S, xPoczatkowe, xKoncowe, xSkok);
	szacowanieBleduAproksymacji_interpolacjaFS3S(daneWezlow_64, liczbaRekordow / 64, nazwyPlikuZapisu.nazwa64_interpolacjaFS3S, daneWezlow, liczbaRekordow);

	//co 16 punkt
	DaneFunkcji *daneWezlow_16 = NULL;
	wyciagnijDaneZPliku_daneOWezlach(&daneWezlow_16, nazwaPlikuDanych, liczbaRekordow, 16);

	iterujPoFunkcji_interpolacjaLagrange(daneWezlow_16, liczbaRekordow / 16, nazwyPlikuZapisu.nazwa16_interpolacjaLagrange, xPoczatkowe, xKoncowe, xSkok);
	szacowanieBleduAproksymacji_interpolacjaLagrange(daneWezlow_16, liczbaRekordow / 16, nazwyPlikuZapisu.nazwa16_interpolacjaLagrange, daneWezlow, liczbaRekordow);

	iterujPoFunkcji_interpolacjaFS3S(daneWezlow_16, liczbaRekordow / 16, nazwyPlikuZapisu.nazwa16_interpolacjaFS3S, xPoczatkowe, xKoncowe, xSkok);
	szacowanieBleduAproksymacji_interpolacjaFS3S(daneWezlow_16, liczbaRekordow / 16, nazwyPlikuZapisu.nazwa16_interpolacjaFS3S, daneWezlow, liczbaRekordow);

	return;
}

int main()
{
	
	//profil WETI-akademik
	NazwyPlikowZapisu nazwyPlikuZapisu_akademik;
	nazwyPlikuZapisu_akademik.nazwa64_interpolacjaLagrange = "../wyniki_akademik_64_Lagrange.txt";
	nazwyPlikuZapisu_akademik.nazwa64_interpolacjaFS3S = "../wyniki_akademik_64_FS3S.txt";
	nazwyPlikuZapisu_akademik.nazwa16_interpolacjaLagrange = "../wyniki_akademik_16_Lagrange.txt";
	nazwyPlikuZapisu_akademik.nazwa16_interpolacjaFS3S = "../wyniki_akademik_16_FS3S.txt";
	wyznaczDaneProfile("../akademik.csv", nazwyPlikuZapisu_akademik, 0.0, 746.7, 0.1);

	//profil Mount Everest
	NazwyPlikowZapisu nazwyPlikuZapisu_mountEverest;
	nazwyPlikuZapisu_mountEverest.nazwa64_interpolacjaLagrange = "../wyniki_mountEverest_64_Lagrange.txt";
	nazwyPlikuZapisu_mountEverest.nazwa64_interpolacjaFS3S = "../wyniki_mountEverest_64_FS3S.txt";
	nazwyPlikuZapisu_mountEverest.nazwa16_interpolacjaLagrange = "../wyniki_mountEverest_16_Lagrange.txt";
	nazwyPlikuZapisu_mountEverest.nazwa16_interpolacjaFS3S = "../wyniki_mountEverest_16_FS3S.txt";
	wyznaczDaneProfile("../mount_everest.csv", nazwyPlikuZapisu_mountEverest, 0.0, 7803.2, 0.5);
	
	//profil Greenwich
	NazwyPlikowZapisu nazwyPlikuZapisu_greenwich;
	nazwyPlikuZapisu_greenwich.nazwa64_interpolacjaLagrange = "../wyniki_greenwich_64_Lagrange.txt";
	nazwyPlikuZapisu_greenwich.nazwa64_interpolacjaFS3S = "../wyniki_greenwich_64_FS3S.txt";
	nazwyPlikuZapisu_greenwich.nazwa16_interpolacjaLagrange = "../wyniki_greenwich_16_Lagrange.txt";
	nazwyPlikuZapisu_greenwich.nazwa16_interpolacjaFS3S = "../wyniki_greenwich_16_FS3S.txt";
	wyznaczDaneProfile("../greenwich.csv", nazwyPlikuZapisu_greenwich, 0.0, 2069.7, 0.1);
	
	//profil kanion Kolorado
	NazwyPlikowZapisu nazwyPlikuZapisu_kanionKolorado;
	nazwyPlikuZapisu_kanionKolorado.nazwa64_interpolacjaLagrange = "../wyniki_kanionKolorado_64_Lagrange.txt";
	nazwyPlikuZapisu_kanionKolorado.nazwa64_interpolacjaFS3S = "../wyniki_kanionKolorado_64_FS3S.txt";
	nazwyPlikuZapisu_kanionKolorado.nazwa16_interpolacjaLagrange = "../wyniki_kanionKolorado_16_Lagrange.txt";
	nazwyPlikuZapisu_kanionKolorado.nazwa16_interpolacjaFS3S = "../wyniki_kanionKolorado_16_FS3S.txt";
	wyznaczDaneProfile("../kanion_kolorado.csv", nazwyPlikuZapisu_kanionKolorado, 0.0, 29870.5, 1.0);

	return 0;
}