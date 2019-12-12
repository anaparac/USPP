Kratke upute o korištenju i pokretanju algoritama.

1. Testiranje stohastickog gradijentnog spusta na podacima o povjerenju korisnika stranice za trgovanje bitcoinom.
- pokrenuti skritpu : "square_hinge" iz komandne linije, nakon što otvorimo bitcoin direktorij unutar matlaba
- algoritam potom ispisuje : 
			- broj unaprijed poznatih podataka
			- broj pogrešnih predznaka
			- relativnu grešku

2. Testiranje stohastickog gradijentnog spusta na podacima o povjerenju korisnika stranice epinions. Podaci : "bitcoin_net.csv".
- pokrenuti skritpu : "square_hinge" iz komandne linije, nakon što otvorimo epinions direktorij unutar matlaba
- algoritam potom ispisuje : 
			- broj unaprijed poznatih podataka
			- broj pogrešnih predznaka
			- relativnu grešku
Upozorenje : zbog velike kolicine podataka izvršavanje traje duže.

Preporuka : Na navedenim tockama 1. i 2. mijenjati parametar k i broj iteracija te uociti promjenu u preciznosti.

3. Testiranje klasteriranja na nasumicno generiranom skupu podataka. Podaci : "epinions_zarez.txt".
- otvoriti datoteku klasteriranje, unutar nje pokrenuti skriptu "Kmeans_test" iz komandne linije.
- algoritam vraca :
		- I - vektor u kojem je a_i element odgovara particiji kojoj i-t element pripada.