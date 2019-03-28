#include <stdio.h>
#include <stdlib.h>

#include <cmath>
#include <fstream>
#include <iostream>
#include <chrono>
#include <opencv2/opencv.hpp>
#include <ctime>
#include <random>

#include "defines.h"
#include "utils.h"

#include "proj_api.h"

#define STEP1_WIN_NAME			  "Heightmap"
#define STEP2_WIN_NAME			  "Edges"
#define ZOOM					  1
#define NUM_MORFOLOGIC_OPERATIONS 2

unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
std::mt19937_64 generator(seed);

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

int bng_1m_accuracy(double x, double y) {
	// e.g. x = -1.8, y = 51.18;  => Easting: 414075.69   Northing: 142326.96

	projPJ pj_src, pj_dst;
	int p;

	// EPSG:4326 definition: http://spatialreference.org/ref/epsg/4326/proj4/
	const char* src = "+proj=krovak +ellps=bessel +towgs84=570.8,85.7,462.8,4.998,1.587,5.261,3.56";
	const char* dst = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs";

	if(!(pj_src = pj_init_plus(src)))
		exit(1);
	if(!(pj_dst = pj_init_plus(dst)))
		exit(1);

	x *= DEG_TO_RAD;
	y *= DEG_TO_RAD;
	p = pj_transform(pj_src, pj_dst, 1, 1, &x, &y, NULL);
	printf("%.2f\t%.2f\n", x, y);
}

/**
 * Mouse clicking callback.
 */
void mouse_probe_handler(int event, int x, int y, int flags, void* param) {
	MouseProbe *probe = (MouseProbe*)param;

	switch(event) {

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
	std::uniform_real_distribution<float> dis(0.0, 255.0);

	return cv::Vec3f(dis(generator), dis(generator), dis(generator));
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
	cv::Mat colorizeMapBuffer = cv::Mat::zeros(edgemap_8uc1_img.size(), edgemap_8uc1_img.type());
	colorizeMapBuffer.at<bool>(x, y) = true;

	std::vector<cv::Point> colorizePixels;
	colorizePixels.push_back(cv::Point(x, y));

	cv::Vec3f randomColor = get_random_color();

	while(!colorizePixels.empty()) {
		cv::Point point = colorizePixels[colorizePixels.size() - 1];

		colorizePixels.pop_back();

		for(int i = 0; i < 4; i++) {
			cv::Point pointToCheck = point + lookForMatrix[i];

			if(pointToCheck.x > 0 && pointToCheck.y > 0 && pointToCheck.x < edgemap_8uc1_img.cols - 1 && pointToCheck.y < edgemap_8uc1_img.rows - 1) {
				bool mapBufferValue = colorizeMapBuffer.at<bool>(pointToCheck.y, pointToCheck.x);

				if(!mapBufferValue && edgemap_8uc1_img.at<uchar>(pointToCheck.y, pointToCheck.x) == edgemap_8uc1_img.at<uchar>(point.y, point.x)) {
					colorizePixels.push_back(pointToCheck);

					colorizeMapBuffer.at<bool>(pointToCheck.y, pointToCheck.x) = true;
				}
			}
		}

		heightmap_show_8uc3_img.at<cv::Vec3b>(point.y, point.x) = randomColor;
	}
}

/**
 * Find the minimum and maximum coordinates in the file.
�* Note that the file is the S-JTSK coordinate system.
 */
void get_min_max(std::vector<const char*> filenames, float *a_min_x, float *a_max_x, float *a_min_y, float *a_max_y, float *a_min_z, float *a_max_z) {
	float x, y, z;
	float min_x, min_y, min_z, max_x, max_y, max_z;
	int l_type;

	min_x = min_y = min_z = INT32_MAX;
	max_x = max_y = max_z = INT32_MIN;

	for(int i = 0; i < filenames.size(); i++) {
		printf("Processing size of image %s\n", filenames.at(i));

		std::ifstream fin(filenames.at(i), std::ios::binary);
		while(fin.read(reinterpret_cast<char*>(&x), sizeof(float))) {
			fin.read(reinterpret_cast<char*>(&y), sizeof(float));
			fin.read(reinterpret_cast<char*>(&z), sizeof(float));

			if(filenames.at(i) == "pt000023.bin") {
				fin.read(reinterpret_cast<char*>(&l_type), sizeof(int));
			}

			min_x = std::min(min_x, x);
			min_y = std::min(min_y, y);
			min_z = std::min(min_z, z);

			max_x = std::max(max_x, x);
			max_y = std::max(max_y, y);
			max_z = std::max(max_z, z);
		}
	}
	
	*a_min_x = min_x;
	*a_min_y = min_y;
	*a_min_z = min_z;

	*a_max_x = max_x;
	*a_max_y = max_y;
	*a_max_z = max_z;

	printf("\n\n");
}


/**
 * Fill the image by data from lidar.
 * All lidar points are stored in a an array that has the dimensions of the image. Then the pixel is assigned
 * a value as an average value range from at the corresponding array element. However, with this simple data access, you will lose data precission.
 * filename - file with binarny data
 * img - input image
 */
void fill_image(std::vector<const char*> filenames, cv::Mat & heightmap_8uc1_img, float min_x, float max_x, float min_y, float max_y, float min_z, float max_z) {
	float fx, fy, fz;
	int l_type;

	// zjistime sirku a vysku obrazu
	int delta_x = round(max_x - min_x + 0.5f);
	int delta_y = round(max_y - min_y + 0.5f);
	int delta_z = round(max_z - min_z + 0.5f);

	heightmap_8uc1_img = cv::Mat::zeros(cv::Size(delta_x + 1, delta_y + 1), CV_8UC1);

	cv::Mat mass = cv::Mat::zeros(cv::Size(delta_x + 1, delta_y + 1), CV_32FC1);
	cv::Mat mass_count = cv::Mat::zeros(cv::Size(delta_x + 1, delta_y + 1), CV_32FC1);

	for(int i = 0; i < filenames.size(); i++) {
		printf("Filling image %s\n", filenames.at(i));

		std::ifstream fin(filenames.at(i), std::ios::binary);

		while(fin.read(reinterpret_cast<char*>(&fx), sizeof(float))) {
			fin.read(reinterpret_cast<char*>(&fy), sizeof(float));
			fin.read(reinterpret_cast<char*>(&fz), sizeof(float));

			if(filenames.at(i) == "pt000023.bin") {
				fin.read(reinterpret_cast<char*>(&l_type), sizeof(int));
			}

			mass.at<float>(heightmap_8uc1_img.rows - (int)round(max_y - fy + 0.5f) - 1, heightmap_8uc1_img.cols - (int)round(max_x - fx + 0.5f) - 1) += fz;
			mass_count.at<int>(heightmap_8uc1_img.rows - (int)round(max_y - fy + 0.5f) - 1, heightmap_8uc1_img.cols - (int)round(max_x - fx + 0.5f) - 1)++;
		}
	}

	mass.forEach<uchar>([&](uchar &pixel, const int *position) {
		int y = position[0], x = position[1];

		float sum = mass.at<float>(y, x);
		float count = mass_count.at<int>(y, x);

		float average = sum / count;

		float fullPercent = max_z - min_z;
		float actualPercent = average - min_z;

		if(count > 0) heightmap_8uc1_img.at<uchar>(y, x) = (actualPercent / fullPercent) * 255.0f;
	});

	printf("\n\n");
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

		if(y > 0 && y < edgemap_8uc1_img.rows &&
			x > 0 && x < edgemap_8uc1_img.cols) {

			if(pixel == 255) {
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
		if(y > 0 && y < edgemap_8uc1_img.rows &&
			x > 0 && x < edgemap_8uc1_img.cols) {

			if(edgemap_8uc1_img.at<uchar>(y - 1, x - 1) != 255 ||
				edgemap_8uc1_img.at<uchar>(y - 1, x + 0) != 255 ||
				edgemap_8uc1_img.at<uchar>(y - 1, x + 1) != 255 ||
				edgemap_8uc1_img.at<uchar>(y + 0, x - 1) != 255 ||
				edgemap_8uc1_img.at<uchar>(y + 0, x + 0) != 255 ||
				edgemap_8uc1_img.at<uchar>(y + 0, x + 1) != 255 ||
				edgemap_8uc1_img.at<uchar>(y + 1, x - 1) != 255 ||
				edgemap_8uc1_img.at<uchar>(y + 1, x + 0) != 255 ||
				edgemap_8uc1_img.at<uchar>(y + 1, x + 1) != 255) {

				tmp_edgemap_8uc1_img.at<uchar>(y, x) = 0;
			} else {
				tmp_edgemap_8uc1_img.at<uchar>(y, x) = 255;
			}
		}
	});

	tmp_edgemap_8uc1_img.copyTo(edgemap_8uc1_img);
}


void write_to_image(cv::Mat input_img, cv::Mat input_height_img, cv::Mat & output_img, const char* file_name) {
	input_img.forEach<cv::Vec3b>([&](cv::Vec3b &pixel, const int *position) {
		int y = position[0], x = position[1];

		uchar height = input_height_img.at<uchar>(y, x);

		output_img.at<cv::Vec4b>(y, x) = cv::Vec4b(pixel.val[0], pixel.val[1], pixel.val[2], (int)height == 0 ? 0 : 255.0f);
	});

	cv::flip(output_img, output_img, 0);

	std::vector<int> compression_params; // Stores the compression parameters

	compression_params.push_back(CV_IMWRITE_PXM_BINARY); // Set to PXM compression
	compression_params.push_back(0); // Set type of PXM in our case PGM

	cv::imwrite(file_name, output_img);
}

void process_lidar(const char *txt_filename, std::vector<const char*> bin_filenames, const char *img_filename) {
	
	MouseProbe *mouse_probe;

	cv::Mat heightmap_8uc1_img;      // image of source of lidar data
	cv::Mat heightmap_show_8uc3_img; // image to detected areas
	cv::Mat edgemap_8uc1_img;        // image for edges

	float absolute_min_x, absolute_min_y, absolute_min_z, absolute_max_x, absolute_max_y, absolute_max_z;
	float absolute_delta_x, absolute_delta_y, absolute_delta_z;
	float min_x, max_x, min_y, max_y, min_z, max_z;
	float delta_x, delta_y, delta_z;

	absolute_min_x = absolute_min_y = absolute_min_z = INT32_MAX;
	absolute_max_x = absolute_max_y = absolute_max_z = INT32_MIN;

	get_min_max(bin_filenames, &min_x, &max_x, &min_y, &max_y, &min_z, &max_z);

	printf("min x: %f, max x: %f\n", min_x, max_x);
	printf("min y: %f, max y: %f\n", min_y, max_y);
	printf("min z: %f, max z: %f\n", min_z, max_z);

	delta_x = max_x - min_x;
	delta_y = max_y - min_y;
	delta_z = max_z - min_z;

	printf("delta x: %f\n", delta_x);
	printf("delta y: %f\n", delta_y);
	printf("delta z: %f\n", delta_z);

	printf("\n\n");

	//bng_1m_accuracy(min_x, min_y);
	//bng_1m_accuracy(max_x, max_y);

	fill_image(bin_filenames, heightmap_8uc1_img, min_x, max_x, min_y, max_y, min_z, max_z);

	cv::resize(heightmap_8uc1_img, heightmap_8uc1_img, cv::Size(heightmap_8uc1_img.cols * 0.2f, heightmap_8uc1_img.rows * 0.2f));

	make_edges(heightmap_8uc1_img, edgemap_8uc1_img);
	binarize_image(edgemap_8uc1_img);

	std::string windowNames[] = { "Show", "Before morfologic operations", "After morfologic operations", "Colorized map" };

	cv::namedWindow(windowNames[0]);
	cv::namedWindow(windowNames[1]);
	cv::namedWindow(windowNames[2]);
	cv::namedWindow(windowNames[3]);

	cv::moveWindow(windowNames[0], 0 * heightmap_8uc1_img.cols, 0 * heightmap_8uc1_img.rows);
	cv::moveWindow(windowNames[1], 1 * edgemap_8uc1_img.cols, 0 * edgemap_8uc1_img.rows);
	cv::moveWindow(windowNames[2], 2 * edgemap_8uc1_img.cols, 0 * edgemap_8uc1_img.rows);
	cv::moveWindow(windowNames[3], 0 * edgemap_8uc1_img.cols, 0 * edgemap_8uc1_img.rows);

	cv::imshow(windowNames[0], heightmap_8uc1_img);
	cv::imshow(windowNames[1], edgemap_8uc1_img);

	auto start = std::chrono::system_clock::now();

	for(int i = 0; i < NUM_MORFOLOGIC_OPERATIONS; i++)
		dilate_and_erode_edgemap(edgemap_8uc1_img);

	auto end = std::chrono::system_clock::now();

	cv::imshow(windowNames[2], edgemap_8uc1_img);

	std::chrono::duration<double> elapsed_seconds = end - start;
	std::time_t end_time = std::chrono::system_clock::to_time_t(end);

	std::cout << "finished computation at " << std::ctime(&end_time)
		<< "elapsed time: " << elapsed_seconds.count() << "s\n";

	mouse_probe = new MouseProbe(heightmap_8uc1_img, heightmap_show_8uc3_img, edgemap_8uc1_img);

	cv::setMouseCallback(windowNames[2], mouse_probe_handler, mouse_probe);

	cv::cvtColor(heightmap_8uc1_img, heightmap_show_8uc3_img, CV_GRAY2BGR);

	//cv::Mat output_img = cv::Mat::zeros(heightmap_8uc1_img.size(), CV_8UC4);

	//write_to_image(heightmap_show_8uc3_img, heightmap_8uc1_img, output_img, "lidar_output.png");

	// wait here for user input using (mouse clicking)
	while(1) {
		cv::imshow(windowNames[3], heightmap_show_8uc3_img);

		int key = cv::waitKey(10);
		if(key == 'q') {
			break;
		}
	}
}


int main(int argc, char *argv[]) {
	std::vector<const char*> bin_files;
	char *txt_file, *bin_file, *img_file;

	if(argc < 4) {
		printf("Not enough command line parameters.\n");
		exit(1);
	}

	bin_files.push_back("pt000001.bin"); // 0
	bin_files.push_back("pt000002.bin"); // 1
	bin_files.push_back("pt000004.bin"); // 2
	bin_files.push_back("pt000005.bin"); // 3
	bin_files.push_back("pt000006.bin"); // 4
	bin_files.push_back("pt000008.bin"); // 5
	bin_files.push_back("pt000009.bin"); // 6
	bin_files.push_back("pt000011.bin"); // 7
	bin_files.push_back("pt000015.bin"); // 8
	bin_files.push_back("pt000016.bin"); // 9
	bin_files.push_back("pt000017.bin"); //10
	bin_files.push_back("pt000018.bin"); //11
	bin_files.push_back("pt000019.bin"); //12
	bin_files.push_back("pt000020.bin"); //13
	bin_files.push_back("pt000023.bin"); //14

	txt_file = argv[1];
	bin_file = (char*)bin_files.at(14);
	img_file = argv[3];

	process_lidar(txt_file, bin_files, img_file);

	delete bin_file;

	return 0;
}
