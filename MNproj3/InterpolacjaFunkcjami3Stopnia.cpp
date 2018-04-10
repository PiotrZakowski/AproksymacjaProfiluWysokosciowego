/*************************************************************
 *	Copyright (c) 2017 Piotr Pawe³ ¯akowski.
 *	All rights reserved.
 *************************************************************/

#include "InterpolacjaFunkcjami3Stopnia.h"

//Interpolacja funkcjami sklejanymi trzeciego stopnia (FS3S).

////Tworzenie macierzy ukladu rownan parametrow M.

/*
 * Zwraca wspolczynnik przekatnej macierzy systemowej.
 */
double lamda(int j, DaneFunkcji *daneWezlow)
{
	double hjplus1 = daneWezlow[j + 1].x - daneWezlow[j].x;
	double hj = daneWezlow[j].x - daneWezlow[j - 1].x;

	return hjplus1 / (hj + hjplus1);
}

/*
 * Zwraca wspolczynnik przekatnej macierzy systemowej.
 */
double micro(int j, DaneFunkcji *daneWezlow)
{
	double hj = daneWezlow[j].x - daneWezlow[j - 1].x;
	double hjplus1 = daneWezlow[j + 1].x - daneWezlow[j].x;

	return hj / (hj + hjplus1);
}

/*
 * Logika tworzenia macierzy systemowej dla interpolacji funkcjami sklejanymi trzeciego stopnia.
 */
void stworzMacierzSystemowa(DaneFunkcji *daneWezlow, double **macierzSystemowa, int N)
{
	//! macierz ma rozmiary (N+1)x(N+1)
	macierzSystemowa[0][0] = 2.0;
	macierzSystemowa[0][1] = 0.0;
	for (int j = 2; j <= N; j++)
		macierzSystemowa[0][j] = 0.0;

	for (int i = 1; i < N; i++)
	{
		for (int j = 0; j <= N; j++)
		{
			if (j == i - 1) //micro
				macierzSystemowa[i][j] = micro(i, daneWezlow);
			else if (j == i)
				macierzSystemowa[i][j] = 2.0;
			else if (j == i + 1) //lambda
				macierzSystemowa[i][j] = lamda(i, daneWezlow);
			else
				macierzSystemowa[i][j] = 0.0;
		}
	}

	for (int j = 0; j <= N - 2; j++)
	{
		macierzSystemowa[N][j] = 0.0;
	}
	macierzSystemowa[N][N - 1] = 0.0;
	macierzSystemowa[N][N] = 2.0;

	return;
}

////Tworzenie wektoru wyrazow wolnych ukladu rownan parametrow M.

/*
 * Zwraca wspolczynnik wektora pobudzenia.
 */
double delta(int j, DaneFunkcji *daneWezlow)
{
	double hj = daneWezlow[j].x - daneWezlow[j - 1].x;
	double hjplus1 = daneWezlow[j + 1].x - daneWezlow[j].x;

	return (6 / (hj + hjplus1))*(((daneWezlow[j + 1].fx - daneWezlow[j].fx) / hjplus1) - ((daneWezlow[j].fx - daneWezlow[j - 1].fx) / hj));
}

/*
 * Logika tworzenia wektora pobudzenia dla interpolacji funkcjami sklejanymi trzeciego stopnia.
 */
void stworzWektorPobudzenia(DaneFunkcji *daneWezlow, double *wektorPobudzenia, int N)
{
	wektorPobudzenia[0] = 0.0;

	for (int j = 1; j <= N - 1; j++)
		wektorPobudzenia[j] = delta(j, daneWezlow);

	wektorPobudzenia[N] = 0.0;

	return;
}

////Rozwiazanie ukladu rownan metoda Gaussa.

/*
 * Tworzy kopie wektora.
 */
void przepiszWektor(double *wektor, double *kopiaWektora, int rozmiar)
{
	for (int i = 0; i < rozmiar; i++)
	{
		kopiaWektora[i] = wektor[i];
	}

	return;
}

/*
 * Tworzy kopie macierzy.
 */
void przepiszMacierz(double **macierz, double **kopiaMacierzy, int rozmiar)
{
	for (int i = 0; i < rozmiar; i++)
	{
		for (int j = 0; j < rozmiar; j++)
		{
			kopiaMacierzy[i][j] = macierz[i][j];
		}
	}

	return;
}

/*
 * Rozwiazanie ukladu rownan liniowych metoda bezposrednia Gaussa.
 */
void algorytmGaussa(double **macierzSystemowa, double *wektorPobudzenia, double *wektorRozwiazan, int rozmiarMacierzy)
{
	double **pomMacierzSystemowa = (double **)malloc(rozmiarMacierzy*sizeof(double **));
	for (int i = 0; i < rozmiarMacierzy; i++)
	{
		pomMacierzSystemowa[i] = (double *)malloc(rozmiarMacierzy*sizeof(double));
	}
	przepiszMacierz(macierzSystemowa, pomMacierzSystemowa, rozmiarMacierzy);

	double *pomWektorPobudzenia = (double *)malloc(rozmiarMacierzy*sizeof(double));
	przepiszWektor(wektorPobudzenia, pomWektorPobudzenia, rozmiarMacierzy);

	//etap postepowania prostego
	double stosunekM = 0.0;
	for (int k = 0; k < rozmiarMacierzy - 1; k++)
	{
		for (int i = k + 1; i < rozmiarMacierzy; i++)
		{
			stosunekM = pomMacierzSystemowa[i][k] / pomMacierzSystemowa[k][k];
			for (int j = 0; j < rozmiarMacierzy; j++)
			{
				pomMacierzSystemowa[i][j] -= stosunekM*pomMacierzSystemowa[k][j];
			}
			pomWektorPobudzenia[i] -= stosunekM*pomWektorPobudzenia[k];
		}
	}

	//etap postepowania odwrotnego
	double suma = 0.0;
	for (int i = rozmiarMacierzy - 1; i >= 0; i--)
	{
		//obliczanie sumy
		suma = 0.0;
		for (int j = i + 1; j < rozmiarMacierzy; j++)
		{
			//suma+=aij*xj
			suma += pomMacierzSystemowa[i][j] * wektorRozwiazan[j];
		}

		//xi=(bi-suma)/aii
		wektorRozwiazan[i] = (pomWektorPobudzenia[i] - suma) / pomMacierzSystemowa[i][i];
	}

	for (int i = 0; i < rozmiarMacierzy; i++)
	{
		free(pomMacierzSystemowa[i]);
	}
	free(pomMacierzSystemowa);
	free(pomWektorPobudzenia);

	return;
}

////Obliczenie wartosci interpolujacej funkcji dla zadanego x.

/*
 * Zwraca wartosci funkcji interpolujacej funkcji sklejanych 3 stopnia w punkcie x.
 */
double interpolacjaFS3S(const DaneFunkcji *daneWezlow, const int liczbaWezlow, const double x, double *wektorParametrowM)
{
	int j;
	for (j = 0; j < liczbaWezlow; j++)
	{
		if (daneWezlow[j].x >= x)
		{
			if (daneWezlow[j].x != 0)
				j--;
			break;
		}
	}

	if (j == liczbaWezlow) //ekstrapolacja
		return 0.0;

	double hjplus1 = daneWezlow[j + 1].x - daneWezlow[j].x;

	double aj = daneWezlow[j].fx;
	double bj = ((daneWezlow[j + 1].fx - daneWezlow[j].fx) / hjplus1) - (((2.0 * wektorParametrowM[j] + wektorParametrowM[j + 1]) / 6.0)*hjplus1);
	double cj = wektorParametrowM[j] / 2.0;
	double dj = (wektorParametrowM[j + 1] - wektorParametrowM[j]) / (6.0*hjplus1);

	return aj + bj*(x - daneWezlow[j].x) + cj*pow(x - daneWezlow[j].x, 2) + dj*pow(x - daneWezlow[j].x, 3);
}

/*
 * Logika obliczania wartosci funkcji interpolujacej funkcji sklejanych 3 stopnia dla kolejnych punktow x nalezacego do <xPoczatkowe;xKoncowej> gdzie x(i+1)=x(i)+xSkok.
 */
void iterujPoFunkcji_interpolacjaFS3S(DaneFunkcji *daneWezlow, int liczbaWezlow, char *nazwaPlikuZapisu, double xPoczatkowe, double xKoncowe, double xSkok)
{
	FILE *plik = fopen(nazwaPlikuZapisu, "w+");

	DaneFunkcji aktualnyPunkt;
	aktualnyPunkt.x = xPoczatkowe;

	int nMacierzy = liczbaWezlow + 1;
	double **macierzSystemowa = (double **)malloc(nMacierzy*sizeof(double **));
	for (int i = 0; i < nMacierzy; i++)
	{
		macierzSystemowa[i] = (double *)malloc(nMacierzy*sizeof(double));
	}
	stworzMacierzSystemowa(daneWezlow, macierzSystemowa, liczbaWezlow);

	double *wektorPobudzenia = (double *)malloc(nMacierzy*sizeof(double));
	stworzWektorPobudzenia(daneWezlow, wektorPobudzenia, liczbaWezlow);

	double *wektorParametrowM = (double *)malloc(nMacierzy*sizeof(double));
	algorytmGaussa(macierzSystemowa, wektorPobudzenia, wektorParametrowM, nMacierzy);

	while (aktualnyPunkt.x < xKoncowe)
	{
		aktualnyPunkt.fx = interpolacjaFS3S(daneWezlow, nMacierzy, aktualnyPunkt.x, wektorParametrowM);

		fprintf(plik, "%lf;%lf\n", aktualnyPunkt.x, aktualnyPunkt.fx);

		aktualnyPunkt.x += xSkok;
	}

	for (int i = 0; i < nMacierzy; i++)
	{
		free(macierzSystemowa[i]);
	}
	free(macierzSystemowa);
	free(wektorPobudzenia);
	free(wektorParametrowM);
	fclose(plik);

	return;
}

////Szacowanie bledu.

/*
 * Szacowanie bledu aproksymacji sredniokwadratowej dyskretnej dla interpolacji funkcji sklejanych 3 stopnia.
 */
void szacowanieBleduAproksymacji_interpolacjaFS3S(DaneFunkcji *daneWezlow, int liczbaWezlow, char *nazwaPlikuZapisu, DaneFunkcji *punktyAproksymowanej, int liczbaPunktowAproksymowanej)
{
	int nMacierzy = liczbaWezlow + 1;
	double **macierzSystemowa = (double **)malloc(nMacierzy*sizeof(double **));
	for (int i = 0; i < nMacierzy; i++)
	{
		macierzSystemowa[i] = (double *)malloc(nMacierzy*sizeof(double));
	}
	stworzMacierzSystemowa(daneWezlow, macierzSystemowa, liczbaWezlow);

	double *wektorPobudzenia = (double *)malloc(nMacierzy*sizeof(double));
	stworzWektorPobudzenia(daneWezlow, wektorPobudzenia, liczbaWezlow);

	double *wektorParametrowM = (double *)malloc(nMacierzy*sizeof(double));
	algorytmGaussa(macierzSystemowa, wektorPobudzenia, wektorParametrowM, nMacierzy);

	FILE *plik = fopen(nazwaPlikuZapisu, "a");

	double suma = 0;
	double interpolowanaWartoscPunktu;
	for (int i = 0; i < liczbaPunktowAproksymowanej; i++)
	{
		interpolowanaWartoscPunktu = interpolacjaFS3S(daneWezlow, nMacierzy, punktyAproksymowanej[i].x, wektorParametrowM);
		suma += pow(punktyAproksymowanej[i].fx - interpolowanaWartoscPunktu, 2);
	}

	fprintf(plik, "Blad aproksymacji interpolacji FS3S:;%lf\n", sqrt(suma));

	for (int i = 0; i < nMacierzy; i++)
	{
		free(macierzSystemowa[i]);
	}
	free(macierzSystemowa);
	free(wektorPobudzenia);
	free(wektorParametrowM);
	fclose(plik);

	return;
}