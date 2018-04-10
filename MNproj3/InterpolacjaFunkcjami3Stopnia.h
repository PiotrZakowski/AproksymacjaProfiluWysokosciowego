/*************************************************************
 *	Copyright (c) 2017 Piotr Pawe³ ¯akowski.
 *	All rights reserved.
 *************************************************************/

#include "Source.h"

//Interpolacja funkcjami sklejanymi trzeciego stopnia (FS3S).

////Obliczenie wartosci interpolujacej funkcji dla zadanego x.

/*
 * Zwraca wartosci funkcji interpolujacej funkcji sklejanych 3 stopnia w punkcie x.
 */
double interpolacjaFS3S(const DaneFunkcji *daneWezlow,const int liczbaWezlow,const double x, double *wektorParametrowM);

/*
 * Logika obliczania wartosci funkcji interpolujacej funkcji sklejanych 3 stopnia dla kolejnych punktow x nalezacego do <xPoczatkowe;xKoncowej> gdzie x(i+1)=x(i)+xSkok.
 */
void iterujPoFunkcji_interpolacjaFS3S(DaneFunkcji *daneWezlow, int liczbaWezlow, char *nazwaPlikuZapisu, double xPoczatkowe, double xKoncowe, double xSkok);

////Szacowanie bledu.

/*
 * Szacowanie bledu aproksymacji sredniokwadratowej dyskretnej dla interpolacji funkcji sklejanych 3 stopnia.
 */
void szacowanieBleduAproksymacji_interpolacjaFS3S(DaneFunkcji *daneWezlow, int liczbaWezlow, char *nazwaPlikuZapisu, DaneFunkcji *punktyAproksymowanej, int liczbaPunktowAproksymowanej);