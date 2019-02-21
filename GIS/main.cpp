#include <stdio.h>
#include <stdlib.h>

#include <cmath>

#include <opencv2/opencv.hpp>

#include "defines.h"
#include "utils.h"

struct MouseProbe {
	cv::Mat & img_;
	cv::Mat & show_img_;
	cv::Mat & edge_img_;

	MouseProbe(cv::Mat & img, cv::Mat & show_img, cv::Mat & edge_img) : img_(img), show_img_(show_img), edge_img_(edge_img)
	{
	}
};

// function declarations
void flood_fill(cv::Mat & src_img, cv::Mat & dst_img, cv::Mat & height_img, const int x, const int y);


/**
 * Callback pro kliknuti mysi.
 * param - obsahuje poiter na libovolna data, ktera mu predame pri vytvareni callbacku
 */
void mouse_probe_handler(int event, int x, int y, int flags, void* param) {
	MouseProbe *probe = (MouseProbe*)param;

	switch (event) {

	case CV_EVENT_LBUTTONDOWN:
		printf("Clicked LEFT at: [ %d, %d ]\n", x, y);
		flood_fill(probe->edge_img_, probe->show_img_, probe->img_, x, y);
		break;

	case CV_EVENT_RBUTTONDOWN:
		printf("Clicked RIGHT at: [ %d, %d ]\n", x, y);
		break;
	}
}


void create_windows(const int width, const int height) {
	cv::namedWindow(STEP1_WIN_NAME, 0);
	cv::namedWindow(STEP2_WIN_NAME, 0);

	cv::resizeWindow(STEP1_WIN_NAME, width*ZOOM, height*ZOOM);
	cv::resizeWindow(STEP2_WIN_NAME, width*ZOOM, height*ZOOM);

} // create_windows


/**
 * Provede flood fill ze zadaneho bodu (x, y) pro okolni body pokud obsahuji stejnou hodnotu,
 * jaka prisla v argumentu value. Funkce rekurzivne vola sama sebe pro sve 4-okoli.
 * src_img - obraz na kterem se bude provadet vyplnovani
 * dst_img - obraz, do ktereho zaznamename vyplneni
 * value - hodnota, pro kterou provedeme vyplneni
 */
void fill_step(cv::Mat & edge_8uc1_img, cv::Mat & heightmap_show_8uc3_img, cv::Mat & heightmap_8uc1_img, const int x, const int y, const uchar value) {
	int width, height;
	int z;

} //fill_step


/**
 * Provede flood fill ze zadaneho bodu (x, y). Funkce si zapamatuje hodnotu na souradnici (x, y)
 * a vyplnuje okoli pomoci funkce fill_step tak dlouho, dokud je hodnota v okolnich bodech stejna.
 * Vyplnovani provadejte na nejakem docasnem obraze, aby nedoslo k poskozeni puvodniho obrazu.
 * src_img - obraz na kterem se bude provadet vyplnovani
 * dst_img - obraz, do ktereho zaznamename vyplneni
 */
void flood_fill(cv::Mat & edge_8uc1_img, cv::Mat & heightmap_show_8uc3_img, cv::Mat & heightmap_8uc1_img, const int x, const int y) {
	cv::Mat tmp_ff_img;

} //flood_fill


/**
 * Zjisti minimalni a maximalni souradnice v zadanem souboru.
 * Nezapomente, ze soubor tvori S-JTSK souradnice.
 */
void get_min_max(const char *filename, float *a_min_x, float *a_max_x, float *a_min_y, float *a_max_y, float *a_min_z, float *a_max_z) {
	FILE *f = NULL;
	float x, y, z;
	float min_x, min_y, min_z, max_x, max_y, max_z;
	int l_type;

	int iteration = 0;

	if ((f = fopen(filename, "r")) == NULL) {
		printf("Read error!");
	}
	else {
		min_x = min_y = min_z = INT32_MAX;
		max_x = max_y = max_z = INT32_MIN;

		while (!feof(f)) {
			fread(&x, sizeof(float), 1, f);
			fread(&y, sizeof(float), 1, f);
			fread(&z, sizeof(float), 1, f);

			fread(&l_type, sizeof(int), 1, f);

			min_x = std::min(min_x, x);
			min_y = std::min(min_y, y);
			min_z = std::min(min_z, z);

			max_x = std::max(max_x, x);
			max_y = std::max(max_y, y);
			max_z = std::max(max_z, z);

			iteration++;
		}
	}

	*a_min_x = min_x;
	*a_min_y = min_y;
	*a_min_z = min_z;

	*a_max_x = max_x;
	*a_max_y = max_y;
	*a_max_z = max_z;

	printf("num of iterations: %d\n", iteration);

	fclose(f);
}


/**
 * Naplni obraz daty z lidaru.
 * Vsechny lidarove body jsou ukladany do pole, ktere ma rozmery obrazu. Pote je jednotlivym pixelum prirazena
 * hodnota jako prumer hodnot z odpovidajiciho prvku pole. Timto jednoduchym pristupem vsak dochazi ke ztrate dat.
 * filename - soubor s binarnimi daty
 * heightmap_8uc1_img - vystupni obrazek
 */
void fill_image(const char *filename, cv::Mat & heightmap_8uc1_img, float min_x, float max_x, float min_y, float max_y, float min_z, float max_z) {
	FILE *f = NULL;
	int delta_x, delta_y, delta_z;
	float fx, fy, fz;
	int x, y, l_type;
	int stride;
	int num_points = 0;
	float range = 0.0f;
	float *sum_height = NULL;
	int *sum_height_count = NULL;

	// zjistime sirku a vysku obrazu
	delta_x = cvRound(max_x - min_x + 0.5f);
	delta_y = cvRound(max_y - min_y + 0.5f);
	delta_z = cvRound(max_z - min_z + 0.5f);

	stride = delta_x;

	// naalokujeme pomocna pole, ve kterych budeme ukaladat hodnoty z lidaru
	// a pocet techto hodnot pro kazdy pixel

	// projdeme soubor a hodnoty priradime do poli

	// hodnoty z pomocneho pole priradime do obrazu
}


void make_edges(const cv::Mat & src_8uc1_img, cv::Mat & edgemap_8uc1_img) {
	cv::Canny(src_8uc1_img, edgemap_8uc1_img, 1, 80);
}


/**
 * Prevede hodnoty obrazu na pouze 2 hodnoty. Hranici mozno nastavit experimentalne.
 */
void binarize_image(cv::Mat & img) {
	int x, y;
	uchar value;

}


/**
 * Provede erozi a dilataci na obraze hran abychom zalepili diry v objektech.
 */
void erode_and_dilate(cv::Mat & edgemap_8uc1_img) {

}


void process_lidar(const char *txt_filename, const char *bin_filename, const char *img_filename) {
	float min_x, max_x, min_y, max_y, min_z, max_z;
	float delta_x, delta_y, delta_z;
	MouseProbe *mouse_probe;

	cv::Mat heightmap_8uc1_img;      // obraz pro vstup lidarovych dat
	cv::Mat heightmap_show_8uc3_img; // obraz pro kresleni nalezenych ploch
	cv::Mat edgemap_8uc1_img;        // obraz pro hrany

	get_min_max(bin_filename, &min_x, &max_x, &min_y, &max_y, &min_z, &max_z);

	printf("min x: %f, max x: %f\n", min_x, max_x);
	printf("min y: %f, max y: %f\n", min_y, max_y);
	printf("min z: %f, max z: %f\n", min_z, max_z);

	delta_x = max_x - min_x;
	delta_y = max_y - min_y;
	delta_z = max_z - min_z;

	printf("delta x: %f\n", delta_x);
	printf("delta y: %f\n", delta_y);
	printf("delta z: %f\n", delta_z);

	// vytvorime obrazky podle informari ze souboru
	/*
	heightmap_8uc1_img = cv::Mat( cvSize( cvRound( delta_x + 0.5f ), cvRound( delta_y + 0.5f ) ), CV_8UC1 );
	heightmap_show_8uc3_img = cv::Mat( cvSize( cvRound( delta_x + 0.5f ), cvRound( delta_y + 0.5f ) ), CV_8UC3 );
	edgemap_8uc1_img = cv::Mat( cvSize( cvRound( delta_x + 0.5f ), cvRound( delta_y + 0.5f ) ), CV_8UC3 );

	create_windows( heightmap_8uc1_img.cols, heightmap_8uc1_img.rows );
	mouse_probe = new MouseProbe( heightmap_8uc1_img, heightmap_show_8uc3_img, edgemap_8uc1_img );

	cv::setMouseCallback( STEP1_WIN_NAME, mouse_probe_handler, mouse_probe );
	cv::setMouseCallback( STEP2_WIN_NAME, mouse_probe_handler, mouse_probe );
	*/

	// naplnime vstupni obraz daty z lidaru
	//fill_image( bin_filename, heightmap_8uc1_img, min_x, max_x, min_y, max_y, min_z, max_z );
	//cv::cvtColor( heightmap_8uc1_img, heightmap_show_8uc3_img, CV_GRAY2RGB );

	// vytvorime obraz hran
	//make_edges( heightmap_8uc1_img, edgemap_8uc1_img );

	// muzeme obraz hran binarizovat, ale v prvni fazi to neni nutne
	//binarize_image( edgemap_8uc1_img );
	//erode_and_dilate( edgemap_8uc1_img );


	// zde cekame na klikani uzivatele
	/*
	cv::imwrite( img_filename, heightmap_8uc1_img );
	while ( 1 ) {
		cv::imshow( STEP1_WIN_NAME, heightmap_show_8uc3_img );
		cv::imshow( STEP2_WIN_NAME, edgemap_8uc1_img );
		int key = cv::waitKey( 10 );
		if ( key == 'q' ) {
			break;
		}
	}
	*/
}


int main(int argc, char *argv[]) {
	char *txt_file, *bin_file, *img_file;

	if (argc < 4) {
		printf("Not enough parameters.\n");
		exit(1);
	}

	txt_file = argv[1];
	bin_file = argv[2];
	img_file = argv[3];

	process_lidar(txt_file, bin_file, img_file);

	system("PAUSE");

	return 0;
}
