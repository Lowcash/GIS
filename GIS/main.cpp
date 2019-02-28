#include <stdio.h>
#include <stdlib.h>

#include <cmath>
#include <fstream>
#include <iostream>

#include <opencv2/opencv.hpp>

#include "defines.h"
#include "utils.h"

#define STEP1_WIN_NAME "Heightmap"
#define STEP2_WIN_NAME "Edges"
#define ZOOM           1


struct MouseProbe {
	cv::Mat & heightmap_8uc1_img_;
	cv::Mat & heightmap_show_8uc3_img_;
	cv::Mat & edgemap_8uc1_img_;

	MouseProbe(cv::Mat & heightmap_8uc1_img, cv::Mat & heightmap_show_8uc3_img, cv::Mat & edgemap_8uc1_img)
		: heightmap_8uc1_img_(heightmap_8uc1_img), heightmap_show_8uc3_img_(heightmap_show_8uc3_img), edgemap_8uc1_img_(edgemap_8uc1_img)
	{
	}
};

// variables

// function declarations
void flood_fill(cv::Mat & src_img, cv::Mat & dst_img, const int x, const int y);


/**
 * Mouse clicking callback.
 */
void mouse_probe_handler(int event, int x, int y, int flags, void* param) {
	MouseProbe *probe = (MouseProbe*)param;

	switch (event) {

	case CV_EVENT_LBUTTONDOWN:
		printf("Clicked LEFT at: [ %d, %d ]\n", x, y);
		flood_fill(probe->edgemap_8uc1_img_, probe->heightmap_show_8uc3_img_, x, y);
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
 * Perform flood fill from the specified point (x, y) for the neighborhood points if they contain the same value,
 * as the one that came in argument 'value'. Function recursicely call itself for its 4-neighborhood.
 *
 * edgemap_8uc1_img - image, in which we perform flood filling
 * heightmap_show_8uc3_img - image, in which we display the filling
 * value - value, for which we'll perform flood filling
 */
void fill_step(cv::Mat & edgemap_8uc1_img, cv::Mat & heightmap_show_8uc3_img, const int x, const int y, const uchar value) {
	int width, height;

} //fill_step


/**
 * Perform flood fill from the specified point (x, y). The function remembers the value at the coordinate (x, y)
 * and fill the neighborhood using 'fill_step' function so long as the value in the neighborhood points are the same.
 * Execute the fill on a temporary image to prevent the original image from being repainted.

 * edgemap_8uc1_img - image, in which we perform flood filling
 * heightmap_show_8uc3_img - image, in which we display the filling
 */
void flood_fill(cv::Mat & edgemap_8uc1_img, cv::Mat & heightmap_show_8uc3_img, const int x, const int y) {
	cv::Mat tmp_edgemap_8uc1_img;

} //flood_fill


/**
 * Find the minimum and maximum coordinates in the file.
 * Note that the file is the S-JTSK coordinate system.
 */
void get_min_max(const char *filename, float *a_min_x, float *a_max_x, float *a_min_y, float *a_max_y, float *a_min_z, float *a_max_z) {
	float x, y, z;
	float min_x, min_y, min_z, max_x, max_y, max_z;
	int l_type;
	
	std::ifstream fin(filename, std::ios::binary);

	min_x = min_y = min_z = INT32_MAX;
	max_x = max_y = max_z = INT32_MIN;

	while (fin.read(reinterpret_cast<char*>(&x), sizeof(float))) {
		fin.read(reinterpret_cast<char*>(&y), sizeof(float));
		fin.read(reinterpret_cast<char*>(&z), sizeof(float));
		fin.read(reinterpret_cast<char*>(&l_type), sizeof(int));

		min_x = std::min(x, min_x);
		min_y = std::min(y, min_y);
		min_z = std::min(z, min_z);

		max_x = std::max(x, max_x);
		max_y = std::max(y, max_y);
		max_z = std::max(z, max_z);
	}

	*a_min_x = min_x;
	*a_min_y = min_y;
	*a_min_z = min_z;

	*a_max_x = max_x;
	*a_max_y = max_y;
	*a_max_z = max_z;
}


/**
 * Fill the image by data from lidar.
 * All lidar points are stored in a an array that has the dimensions of the image. Then the pixel is assigned
 * a value as an average value range from at the corresponding array element. However, with this simple data access, you will lose data precission.
 * filename - file with binarny data
 * img - input image
 */
void fill_image(const char *filename, cv::Mat & heightmap_8uc1_img, float min_x, float max_x, float min_y, float max_y, float min_z, float max_z) {
	FILE *f = NULL;
	int delta_x, delta_y, delta_z;
	float fx, fy, fz;
	int l_type;
	int stride;
	int num_points = 1;
	float range = 0.0f;
	float *sum_height = NULL;
	int *sum_height_count = NULL;

	// zjistime sirku a vysku obrazu
	delta_x = round(max_x - min_x + 0.5f);
	delta_y = round(max_y - min_y + 0.5f);
	delta_z = round(max_z - min_z + 0.5f);

	stride = delta_x;

	heightmap_8uc1_img = cv::Mat(cv::Size(delta_x + 1, delta_y + 1), CV_8UC1);

	cv::Mat mass_count;
	heightmap_8uc1_img.copyTo(mass_count);

	std::ifstream fin(filename, std::ios::binary);

	// 1:
	// We allocate helper arrays, in which we store values from the lidar
	// and the number of these values for each pixel

	// 2:
	// go through the file and assign values to the field
	// beware that in S-JTSK the beginning of the co-ordinate system is at the bottom left,
	// while in the picture it is top left
	while (fin.read(reinterpret_cast<char*>(&fx), sizeof(float))) {
		fin.read(reinterpret_cast<char*>(&fy), sizeof(float));
		fin.read(reinterpret_cast<char*>(&fz), sizeof(float));
		fin.read(reinterpret_cast<char*>(&l_type), sizeof(int));

		heightmap_8uc1_img.at<uchar>(heightmap_8uc1_img.rows - (int)round(max_y - fy + 0.5f) - 1, heightmap_8uc1_img.cols - (int)round(max_x - fx + 0.5f) - 1) += (int)round(max_z - fz + 0.5f) - 1;
		mass_count.at<uchar>(heightmap_8uc1_img.rows - (int)round(max_y - fy + 0.5f) - 1, heightmap_8uc1_img.cols - (int)round(max_x - fx + 0.5f) - 1)++;
	}
	
	for (int y = 0; y < heightmap_8uc1_img.rows; y++) {
		for (int x = 0; x < heightmap_8uc1_img.cols; x++) {
			float count = mass_count.at<uchar>(y, x);
			
			if (count > 0) {
				//heightmap_8uc1_img.at<uchar>(y, x) /= count;

				heightmap_8uc1_img.at<uchar>(y, x) = heightmap_8uc1_img.at<uchar>(y, x) / max_z * 255.0f;
			}
		}
	}
	// 3:
	// assign values from the helper field into the image

}


void make_edges(const cv::Mat & src_8uc1_img, cv::Mat & edgemap_8uc1_img) {
	cv::Canny(src_8uc1_img, edgemap_8uc1_img, 1, 80);
}


/**
 * Transforms the image so it contains only two values.
 * Threshold may be set experimentally.
 */
void binarize_image(cv::Mat & src_8uc1_img) {
	int x, y;
	uchar value;

}


void dilate_and_erode_edgemap(cv::Mat & edgemap_8uc1_img) {
}


void process_lidar(const char *txt_filename, const char *bin_filename, const char *img_filename) {
	float min_x, max_x, min_y, max_y, min_z, max_z;
	float delta_x, delta_y, delta_z;
	MouseProbe *mouse_probe;

	cv::Mat heightmap_8uc1_img;      // image of source of lidar data
	cv::Mat heightmap_show_8uc3_img; // image to detected areas
	cv::Mat edgemap_8uc1_img;        // image for edges

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

	fill_image(bin_filename, heightmap_8uc1_img, min_x, max_x, min_y, max_y, min_z, max_z);

	cv::imshow("Result", heightmap_8uc1_img);
	cv::waitKey(0);

	// create images according to data from the source file
	/*
	heightmap_8uc1_img = cv::Mat( cvSize( cvRound( delta_x + 0.5f ), cvRound( delta_y + 0.5f ) ), CV_8UC1 );
	heightmap_show_8uc3_img = cv::Mat( cvSize( cvRound( delta_x + 0.5f ), cvRound( delta_y + 0.5f ) ), CV_8UC3 );
	edgemap_8uc1_img = cv::Mat( cvSize( cvRound( delta_x + 0.5f ), cvRound( delta_y + 0.5f ) ), CV_8UC3 );

	create_windows( heightmap_8uc1_img.cols, heightmap_8uc1_img.rows );
	mouse_probe = new MouseProbe( heightmap_8uc1_img, heightmap_show_8uc3_img, edgemap_8uc1_img );

	cv::setMouseCallback( STEP1_WIN_NAME, mouse_probe_handler, mouse_probe );
	cv::setMouseCallback( STEP2_WIN_NAME, mouse_probe_handler, mouse_probe );

	printf( "Image w=%d, h=%d\n", heightmap_8uc1_img.cols, heightmap_8uc1_img.rows );
	*/

	// fill the image with data from lidar scanning
	//fill_image( bin_filename, heightmap_8uc1_img, min_x, max_x, min_y, max_y, min_z, max_z );
	//cv::cvtColor( heightmap_8uc1_img, heightmap_show_8uc3_img, CV_GRAY2RGB );

	// create edge map from the height image
	//make_edges( heightmap_8uc1_img, edgemap_8uc1_img );

	// binarize image, so we can easily process it in the next step
	//binarize_image( edgemap_8uc1_img );

	// implement image dilatation and erosion
	//dilate_and_erode_edgemap( edgemap_8uc1_img );

	//cv::imwrite( img_filename, heightmap_8uc1_img );

	// wait here for user input using (mouse clicking)
	/*
	while ( 1 ) {
		cv::imshow( STEP1_WIN_NAME, heightmap_show_8uc3_img );
		//cv::imshow( STEP2_WIN_NAME, edgemap_8uc1_img );
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
		printf("Not enough command line parameters.\n");
		exit(1);
	}

	txt_file = argv[1];
	bin_file = argv[2];
	img_file = argv[3];

	process_lidar(txt_file, bin_file, img_file);

	return 0;
}
