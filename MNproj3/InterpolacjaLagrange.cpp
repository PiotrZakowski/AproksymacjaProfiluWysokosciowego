/*************************************************************
 *	Copyright (c) 2017 Piotr Pawe³ ¯akowski.
 *	All rights reserved.
 *************************************************************/

#include "InterpolacjaLagrange.h"

//Interpolacja Lagrange'a.

////Obliczenie wartosci interpolujacej funkcji dla zadanego x.

/*
 * Zwraca wartosci funkcji interpolujacej Lagrange'a w punkcie x.
 */
double interpolacjaLagrange(const DaneFunkcji *daneWezlow, const int liczbaWezlow, const double x )
{
	double wynik = 0;
	double P, PLicznik, PMianownik;
	for (int k = 0; k < liczbaWezlow; k++)
	{
		PLicznik = 1;
		for (int i = 0; i < liczbaWezlow; i++)
		{
			if (k != i)
				PLicznik *= x - daneWezlow[i].x;
		}

		PMianownik = 1;
		for (int j = 0; j < liczbaWezlow; j++)
		{
			if (k != j)
				PMianownik *= daneWezlow[k].x - daneWezlow[j].x;
		}

		P = PLicznik / PMianownik;
		wynik += daneWezlow[k].fx*P;
	}

	return wynik;
}

/*
 * Logika obliczania wartosci funkcji interpolujacej Lagrange'a dla kolejnych punktow x nalezacego do <xPoczatkowe;xKoncowej> gdzie x(i+1)=x(i)+xSkok.
 */
void iterujPoFunkcji_interpolacjaLagrange(DaneFunkcji *daneWezlow, int liczbaWezlow, char *nazwaPlikuZapisu, double xPoczatkowe, double xKoncowe, double xSkok)
{
	FILE *plik = fopen(nazwaPlikuZapisu, "w+");

	DaneFunkcji aktualnyPunkt;
	aktualnyPunkt.x = xPoczatkowe;

	while (aktualnyPunkt.x < xKoncowe)
	{
		aktualnyPunkt.fx = interpolacjaLagrange(daneWezlow, liczbaWezlow, aktualnyPunkt.x);

		fprintf(plik, "%lf;%lf\n", aktualnyPunkt.x, aktualnyPunkt.fx);

		aktualnyPunkt.x += xSkok;
	}

	fclose(plik);

	return;
}

/*
 * Szacowanie bledu aproksymacji sredniokwadratowej dyskretnej dla interpolacji Lagrange'a.
 */
void szacowanieBleduAproksymacji_interpolacjaLagrange(DaneFunkcji *daneWezlow, int liczbaWezlow, char *nazwaPlikuZapisu, DaneFunkcji *punktyAproksymowanej, int liczbaPunktowAproksymowanej)
{
	FILE *plik = fopen(nazwaPlikuZapisu, "a");

	double suma = 0;
	double interpolowanaWartoscPunktu;
	for (int i = 0; i < liczbaPunktowAproksymowanej; i++)
	{
		interpolowanaWartoscPunktu = interpolacjaLagrange(daneWezlow, liczbaWezlow, punktyAproksymowanej[i].x);
		suma += pow(punktyAproksymowanej[i].fx - interpolowanaWartoscPunktu, 2);
	}

	fprintf(plik, "Blad aproksymacji interpolacji Lagrange'a:;%lf\n", sqrt(suma));
	
	fclose(plik);

	return;
}