#include <stdio.h>
#include <stdlib.h>

#include <cmath>
#include <fstream>
#include <iostream>
#include <chrono>
#include <opencv2/opencv.hpp>
#include <ctime>

#include "defines.h"
#include "utils.h"

#define STEP1_WIN_NAME			  "Heightmap"
#define STEP2_WIN_NAME			  "Edges"
#define ZOOM					  1
#define NUM_MORFOLOGIC_OPERATIONS 2

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

cv::Vec3f get_random_color() {
	return cv::Vec3f((rand() % 256 + 1), (rand() % 256 + 1), (rand() % 256 + 1));
}

cv::Point lookForMatrix[]{
						cv::Point(1, 0),
		cv::Point(0, 1),				cv::Point(0, -1),
						cv::Point(-1, 0)
};

/**
 * Perform flood fill from the specified point (x, y). The function remembers the value at the coordinate (x, y)
 * and fill the neighborhood using 'fill_step' function so long as the value in the neighborhood points are the same.
 * Execute the fill on a temporary image to prevent the original image from being repainted.

 * edgemap_8uc1_img - image, in which we perform flood filling
 * heightmap_show_8uc3_img - image, in which we display the filling
 */
void flood_fill(cv::Mat & edgemap_8uc1_img, cv::Mat & heightmap_show_8uc3_img, const int x, const int y) {
	const int reverse_x = y;
	const int reverse_y = x;

	cv::Mat colorizeMapBuffer = cv::Mat::zeros(edgemap_8uc1_img.size(), edgemap_8uc1_img.type());
	colorizeMapBuffer.at<uchar>(reverse_y, reverse_x) = 255;

	std::vector<cv::Point> colorizePixels;
	colorizePixels.push_back(cv::Point(reverse_y, reverse_x));
	
	cv::Vec3f randomColor = get_random_color();

	while (!colorizePixels.empty()) {
		cv::Point point = colorizePixels[colorizePixels.size() - 1];
		
		colorizePixels.pop_back();

		for (int i = 0; i < 4; i++) {
			cv::Point pointToCheck = point + lookForMatrix[i];

			if (pointToCheck.x > 0 && pointToCheck.y > 0 && pointToCheck.x < edgemap_8uc1_img.cols - 1 && pointToCheck.y < edgemap_8uc1_img.rows - 1) {
				uchar mapBufferValue = colorizeMapBuffer.at<uchar>(pointToCheck.y, pointToCheck.x);

				if (mapBufferValue != 255 && edgemap_8uc1_img.at<uchar>(pointToCheck.y, pointToCheck.x) == edgemap_8uc1_img.at<uchar>(point.y, point.x)) {
					colorizePixels.push_back(pointToCheck);

					colorizeMapBuffer.at<uchar>(pointToCheck.y, pointToCheck.x) = 255;
				}
			}
		}

		heightmap_show_8uc3_img.at<cv::Vec3b>(point.y, point.x) = randomColor;
	}
}


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
	float fx, fy, fz;
	int l_type;

	// zjistime sirku a vysku obrazu
	int delta_x = round(max_x - min_x + 0.5f);
	int delta_y = round(max_y - min_y + 0.5f);
	int delta_z = round(max_z - min_z + 0.5f);

	heightmap_8uc1_img = cv::Mat::zeros(cv::Size(delta_x + 1, delta_y + 1), CV_8UC1);

	cv::Mat mass = cv::Mat::zeros(cv::Size(delta_x + 1, delta_y + 1), CV_32FC1);
	cv::Mat mass_count = cv::Mat::zeros(cv::Size(delta_x + 1, delta_y + 1), CV_32FC1);

	std::ifstream fin(filename, std::ios::binary);

	while (fin.read(reinterpret_cast<char*>(&fx), sizeof(float))) {
		fin.read(reinterpret_cast<char*>(&fy), sizeof(float));
		fin.read(reinterpret_cast<char*>(&fz), sizeof(float));
		fin.read(reinterpret_cast<char*>(&l_type), sizeof(int));

		mass.at<float>(heightmap_8uc1_img.rows - (int)round(max_y - fy + 0.5f) - 1, heightmap_8uc1_img.cols - (int)round(max_x - fx + 0.5f) - 1) += fz;
		mass_count.at<int>(heightmap_8uc1_img.rows - (int)round(max_y - fy + 0.5f) - 1, heightmap_8uc1_img.cols - (int)round(max_x - fx + 0.5f) - 1)++;
	}

	mass.forEach<uchar>([&](uchar &pixel, const int *position) {
		int y = position[0], x = position[1];

		float sum = mass.at<float>(y, x);
		float count = mass_count.at<int>(y, x);

		float average = sum / count;

		float fullPercent = max_z - min_z;
		float actualPercent = average - min_z;

		if (count > 0) {
			heightmap_8uc1_img.at<uchar>(y, x) = (actualPercent / fullPercent) * 255.0f;
		}
	});
}


void make_edges(const cv::Mat & src_8uc1_img, cv::Mat & edgemap_8uc1_img) {
	cv::Canny(src_8uc1_img, edgemap_8uc1_img, 1, 80);
}


/**
 * Transforms the image so it contains only two values.
 * Threshold may be set experimentally.
 */
void binarize_image(cv::Mat & src_8uc1_img) {
	src_8uc1_img.forEach<unsigned char>([&](uchar &pixel, const int *position) {
		pixel = pixel > 127 ? 255 : 0;
	});
}


void dilate_and_erode_edgemap(cv::Mat & edgemap_8uc1_img) {
	cv::Mat tmp_edgemap_8uc1_img = cv::Mat::zeros(edgemap_8uc1_img.size(), edgemap_8uc1_img.type());

	// dilate
	edgemap_8uc1_img.forEach<uchar>([&](uchar &pixel, const int *position) {
		int y = position[0], x = position[1];

		if (y > 0 && y < edgemap_8uc1_img.rows &&
			x > 0 && x < edgemap_8uc1_img.cols) {

			if (pixel == 255) {
				tmp_edgemap_8uc1_img.at<uchar>(y - 1, x - 1) = 255;
				tmp_edgemap_8uc1_img.at<uchar>(y - 1, x + 0) = 255;
				tmp_edgemap_8uc1_img.at<uchar>(y - 1, x + 1) = 255;
				tmp_edgemap_8uc1_img.at<uchar>(y + 0, x - 1) = 255;
				tmp_edgemap_8uc1_img.at<uchar>(y + 0, x + 0) = 255;
				tmp_edgemap_8uc1_img.at<uchar>(y + 0, x + 1) = 255;
				tmp_edgemap_8uc1_img.at<uchar>(y + 1, x - 1) = 255;
				tmp_edgemap_8uc1_img.at<uchar>(y + 1, x + 0) = 255;
				tmp_edgemap_8uc1_img.at<uchar>(y + 1, x + 1) = 255;
			}
		}
	});

	tmp_edgemap_8uc1_img.copyTo(edgemap_8uc1_img);

	tmp_edgemap_8uc1_img.zeros(edgemap_8uc1_img.size(), edgemap_8uc1_img.type());

	// erode
	edgemap_8uc1_img.forEach<unsigned char>([&](uchar &pixel, const int *position) {
		int y = position[0], x = position[1];
		if (y > 0 && y < edgemap_8uc1_img.rows &&
			x > 0 && x < edgemap_8uc1_img.cols) {

			if (edgemap_8uc1_img.at<uchar>(y - 1, x - 1) != 255 ||
				edgemap_8uc1_img.at<uchar>(y - 1, x + 0) != 255 ||
				edgemap_8uc1_img.at<uchar>(y - 1, x + 1) != 255 ||
				edgemap_8uc1_img.at<uchar>(y + 0, x - 1) != 255 ||
				edgemap_8uc1_img.at<uchar>(y + 0, x + 0) != 255 ||
				edgemap_8uc1_img.at<uchar>(y + 0, x + 1) != 255 ||
				edgemap_8uc1_img.at<uchar>(y + 1, x - 1) != 255 ||
				edgemap_8uc1_img.at<uchar>(y + 1, x + 0) != 255 ||
				edgemap_8uc1_img.at<uchar>(y + 1, x + 1) != 255) {

				tmp_edgemap_8uc1_img.at<uchar>(y, x) = 0;
			}
			else {
				tmp_edgemap_8uc1_img.at<uchar>(y, x) = 255;
			}
		}
	});

	tmp_edgemap_8uc1_img.copyTo(edgemap_8uc1_img);
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

	make_edges(heightmap_8uc1_img, edgemap_8uc1_img);
	binarize_image(edgemap_8uc1_img);
	
	std::string windowNames[] = { "Show", "Before morfologic operations", "After morfologic operations", "Colorized map" };

	cv::namedWindow(windowNames[0]);
	cv::namedWindow(windowNames[1]);
	cv::namedWindow(windowNames[2]);
	cv::namedWindow(windowNames[3]);

	cv::moveWindow(windowNames[0], 0 * heightmap_8uc1_img.cols, 0 * heightmap_8uc1_img.rows);
	cv::moveWindow(windowNames[1], 1 * edgemap_8uc1_img.cols  , 0 * edgemap_8uc1_img.rows);
	cv::moveWindow(windowNames[2], 2 * edgemap_8uc1_img.cols  , 0 * edgemap_8uc1_img.rows);
	cv::moveWindow(windowNames[3], 0 * edgemap_8uc1_img.cols, 1 * edgemap_8uc1_img.rows);

	cv::imshow(windowNames[0], heightmap_8uc1_img);
	cv::imshow(windowNames[1], edgemap_8uc1_img);

	auto start = std::chrono::system_clock::now();

	for (int i = 0; i < NUM_MORFOLOGIC_OPERATIONS; i++) {
		dilate_and_erode_edgemap(edgemap_8uc1_img);
	}

	auto end = std::chrono::system_clock::now();

	std::chrono::duration<double> elapsed_seconds = end - start;
	std::time_t end_time = std::chrono::system_clock::to_time_t(end);

	std::cout << "finished computation at " << std::ctime(&end_time)
		<< "elapsed time: " << elapsed_seconds.count() << "s\n";

	cv::imshow(windowNames[2], edgemap_8uc1_img);

	heightmap_show_8uc3_img = cv::Mat::zeros(heightmap_8uc1_img.size(), CV_8UC3);

	heightmap_8uc1_img.forEach<unsigned char>([&](uchar &pixel, const int *position) {
		int y = position[0], x = position[1];

		heightmap_show_8uc3_img.at<cv::Vec3b>(y, x) = cv::Vec3b(pixel, pixel, pixel);
	});

	mouse_probe = new MouseProbe(heightmap_8uc1_img, heightmap_show_8uc3_img, edgemap_8uc1_img);

	cv::setMouseCallback(windowNames[2], mouse_probe_handler, mouse_probe);

	// wait here for user input using (mouse clicking)
	while ( 1 ) {
		cv::imshow(windowNames[3], heightmap_show_8uc3_img );
		
		int key = cv::waitKey( 10 );
		if ( key == 'q' ) {
			break;
		}
	}
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
