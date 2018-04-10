/*************************************************************
 *	Copyright (c) 2017 Piotr Pawe³ ¯akowski.
 *	All rights reserved.
 *************************************************************/

#include "Source.h"

//Interpolacja Lagrange'a.

////Obliczenie wartosci interpolujacej funkcji dla zadanego x.

/*
 * Zwraca wartosci funkcji interpolujacej Lagrange'a w punkcie x.
 */
double interpolacjaLagrange(const DaneFunkcji *daneWezlow, const int liczbaWezlow, const double x);

/*
 * Logika obliczania wartosci funkcji interpolujacej Lagrange'a dla kolejnych punktow x nalezacego do <xPoczatkowe;xKoncowej> gdzie x(i+1)=x(i)+xSkok.
 */
void iterujPoFunkcji_interpolacjaLagrange(DaneFunkcji *daneWezlow, int liczbaWezlow, char *nazwaPlikuZapisu, double xPoczatkowe, double xKoncowe, double xSkok);

////Szacowanie bledu.

/*
 * Szacowanie bledu aproksymacji sredniokwadratowej dyskretnej dla interpolacji Lagrange'a.
 */
void szacowanieBleduAproksymacji_interpolacjaLagrange(DaneFunkcji *daneWezlow, int liczbaWezlow, char *nazwaPlikuZapisu, DaneFunkcji *punktyAproksymowanej, int liczbaPunktowAproksymowanej);