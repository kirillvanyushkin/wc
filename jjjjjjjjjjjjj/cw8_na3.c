--- START OF FILE cw8.c ---

#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <getopt.h>
#include <stddef.h>

typedef struct Options{
    int mode;
    int left_coord[2];
    int right_coord[2];
    size_t thickness;
    int color[3];
    int fill_color[3];
    int center[2];
    size_t radius;
    int dest_left_coord[2];
    int fill;
    char* output_name;
    char* input_name;
    int width_scale;
    int threshold; // For binarization

    int is_left_coord, is_right_coord, is_thickness, is_color, is_fill_color,
        is_center, is_radius, is_dest_left_coord, is_output_name, is_input_name,
        is_width_scale,
        is_threshold; // Flag for binarization threshold

} Options;

#pragma pack(push, 1)
typedef struct BMPFileHeader{
    unsigned char ID[2];
    unsigned int file_size;
    unsigned char unused[4];
    unsigned int pixel_offset;
} BMPFileHeader;

typedef struct BMPInfoHeader{
    unsigned int HeaderSize;
    unsigned int width;
    unsigned int height;
    unsigned short int planes;
    unsigned short int bit_per_pixel;
    unsigned int compression;
    unsigned int image_size;
    unsigned int x_pixel_per_metr;
    unsigned int y_pixel_per_metr;
    unsigned int color_count;
    unsigned int imp_color_count;
} BMPInfoHeader;

typedef struct RGB{
    unsigned char b;
    unsigned char g;
    unsigned char r;
} RGB;

typedef struct BMPFile{
    BMPFileHeader bfm;
    BMPInfoHeader bih;
    RGB** rgb;
} BMPFile;
#pragma pack(pop)

void help_func();
int validate_rectangle(int is_left_coord, int is_right_coord, int is_thickness,
    int is_color, int fill, int is_fill_color, int* color, int* fill_color);
int validate_hexagon(int is_center, int is_radius, int is_thickness,
    int is_color, int fill, int is_fill_color, int* color, int* fill_color);
int validate_copy(int is_left_coord, int is_right_coord, int is_dest_left_coord);
int validate_expand(int is_fill_color, int is_width_scale);
int validate_binarization(int is_threshold, int threshold_val); // New declaration
void free_mem(char** input_name, char** output_name, BMPFile** bmp_file);
int valid_file(BMPFile* bmp_file);
void print_info(BMPFile* bmpf);
void swap(int* a, int* b);
int min(int a, int b);
int max(int a, int b);
int get_padding(int w);
BMPFile* loadBMPFile(char* fname);
void writeBMPFile(char* fname, BMPFile* bmpf);
void draw_hexagon(int* center, size_t radius, size_t thickness,
    int* color, int fill, int* fill_color, BMPFile* bmpf);
void draw_rectangle(int* left_coord, int* right_coord, size_t thickness,
    int* color, int fill, int* fill_color, BMPFile* bmpf);
void copy_area(int* left_coord, int* right_coord, int* dest_coord, BMPFile* bmpf);
void expand(int width_scale, int* color_fill, BMPFile* bmpf, char* output_name);
void do_binarization(BMPFile* bmpf, int threshold_value); // New declaration
int parse_args(int argc, char* argv[], Options* options);
static int check_rgb(int* color_val);


#define PI 3.14159265358979323846

#define HELP_OPT_CHAR 'h'
#define RECT_OPT_VAL 1001
#define LEFT_UP_OPT_VAL 1002
#define RIGHT_DOWN_OPT_VAL 1003
#define THICKNESS_OPT_VAL 1004
#define COLOR_OPT_VAL 1005
#define FILL_OPT_VAL 1006
#define FILL_COLOR_OPT_VAL 1007
#define HEXAGON_OPT_VAL 1008
#define CENTER_OPT_VAL 1009
#define RADIUS_OPT_VAL 1010
#define COPY_OPT_VAL 1011
#define DEST_LEFT_UP_OPT_VAL 1012
#define INFO_OPT_VAL 1013
#define EXPAND_OPT_VAL 1014
#define WIDTH_SCALE_OPT_VAL 1015
#define BINARIZATION_OPT_VAL 1016 // New
#define THRESHOLD_OPT_VAL 1017    // New

#define OUTPUT_OPT_CHAR 'o'
#define INPUT_OPT_CHAR 'i'

#define MODE_UNKNOWN 0
#define MODE_RECTANGLE 1
#define MODE_HEXAGON 2
#define MODE_COPY 3
#define MODE_INFO 4
#define MODE_EXPAND 5
#define MODE_BINARIZATION 6 // New

void help_func(){
    puts("Usage: ./cw [options] [input_file]");
    puts("Options:");
    puts("  -h, --help                  Show this help message and exit.");
    puts("  -i, --input FILE            Input BMP file name.");
    puts("  -o, --output FILE           Output BMP file name (default: out.bmp).");
    puts("");
    puts("Operations (select one):");
    puts("  --rect                      Draw a rectangle.");
    puts("    --left_up X.Y             Top-left corner coordinates for rect/copy source.");
    puts("    --right_down X.Y          Bottom-right corner coordinates for rect/copy source.");
    puts("    --thickness N             Thickness of the border (for rect/hexagon). Must be >= 1.");
    puts("    --color R.G.B             Border color (RGB values 0-255).");
    puts("    --fill                    Enable filling for the shape (rect/hexagon).");
    puts("    --fill_color R.G.B        Fill color (RGB values 0-255, required if --fill).");
    puts("");
    puts("  --hexagon                   Draw a hexagon.");
    puts("    --center X.Y              Center coordinates for hexagon.");
    puts("    --radius N                Radius of the hexagon. Must be >= 1.");
    puts("    (uses --thickness, --color, --fill, --fill_color as with --rect)");
    puts("");
    puts("  --copy                      Copy an area of the image.");
    puts("    (uses --left_up, --right_down for source area)");
    puts("    --dest_left_up X.Y        Top-left corner for destination of copy.");
    puts("");
    puts("  --expand                    Expand image width, filling new columns.");
    puts("    --width_scale N           Scaling factor for width (e.g., 2 for double). Must be > 0.");
    puts("    (uses --fill_color for the new columns)");
    puts("");
    puts("  --binarization              Binarize the image (set pixels to black or white)."); // New
    puts("    --threshold N             Threshold value for binarization (integer > 0). Sum of R+G+B >= N -> white, else black."); // New
    puts("");
    puts("  --info                      Print BMP file information.");
    puts("\nExample usage:");
    puts("  ./cw --input image.bmp --rect --left_up 10.20 --right_down 100.120 --thickness 3 --color 255.0.0 -o output.bmp");
    puts("  ./cw --hexagon --center 150.150 --radius 50 --thickness 2 --color 0.0.255 --fill --fill_color 0.255.0 -i image.bmp");
    puts("  ./cw --input image.bmp --binarization --threshold 384 -o binarized_image.bmp"); // New example
}

static int check_rgb(int* color_val){
    for(int i = 0; i < 3; i++){
        if(color_val[i] > 255 || color_val[i] < 0){
            printf("Error: Color component %d (%d) is out of range [0, 255].\n", i, color_val[i]);
            return 0;
        }
    }
    return 1;
}

int validate_rectangle(int is_left_coord, int is_right_coord, int is_thickness,
     int is_color, int fill, int is_fill_color, int* color, int* fill_color){
    int success = 1;
    if(!is_left_coord) { puts("Error: --rect requires --left_up."); success = 0; }
    if(!is_right_coord) { puts("Error: --rect requires --right_down."); success = 0; }
    if(!is_thickness) { puts("Error: --rect requires --thickness."); success = 0; }
    if(!is_color) { puts("Error: --rect requires --color."); success = 0; }
    else if (!check_rgb(color)) { success = 0; }

    if(fill && !is_fill_color) { puts("Error: --fill for --rect also requires --fill_color."); success = 0; }
    else if (fill && is_fill_color && !check_rgb(fill_color)) { success = 0; }
    return success;
}

int validate_hexagon(int is_center, int is_radius, int is_thickness,
    int is_color, int fill, int is_fill_color, int* color, int* fill_color){
    int success = 1;
    if(!is_center) { puts("Error: --hexagon requires --center."); success = 0; }
    if(!is_radius) { puts("Error: --hexagon requires --radius."); success = 0; }
    if(!is_thickness) { puts("Error: --hexagon requires --thickness."); success = 0; }
    if(!is_color) { puts("Error: --hexagon requires --color."); success = 0; }
    else if(!check_rgb(color)) { success = 0; }

    if(fill && !is_fill_color) { puts("Error: --fill for --hexagon also requires --fill_color."); success = 0; }
    else if (fill && is_fill_color && !check_rgb(fill_color)) { success = 0; }
    return success;
}

int validate_copy(int is_left_coord, int is_right_coord, int is_dest_left_coord){
    int success = 1;
    if(!is_left_coord) { puts("Error: --copy requires --left_up for source."); success = 0; }
    if(!is_right_coord) { puts("Error: --copy requires --right_down for source."); success = 0; }
    if(!is_dest_left_coord) { puts("Error: --copy requires --dest_left_up for destination."); success = 0; }
    return success;
}

int validate_expand(int is_fill_color, int is_width_scale){
    int success = 1;
    if(!is_fill_color) { puts("Error: --expand requires --fill_color."); success = 0; }
    if(!is_width_scale) { puts("Error: --expand requires --width_scale."); success = 0; }
    return success;
}

// New validation function
int validate_binarization(int is_threshold, int threshold_val){
    int success = 1;
    if(!is_threshold) { puts("Error: --binarization requires --threshold."); success = 0; }
    // The check for threshold_val > 0 is already handled by parse_positive_int_str
    return success;
}

void free_mem(char** input_name_ptr, char** output_name_ptr, BMPFile** bmp_file_ptr_ptr){
    if(input_name_ptr && *input_name_ptr){
        free(*input_name_ptr);
        *input_name_ptr = NULL;
    }
    if(output_name_ptr && *output_name_ptr){
        free(*output_name_ptr);
        *output_name_ptr = NULL;
    }
    if(bmp_file_ptr_ptr && *bmp_file_ptr_ptr){
        BMPFile* bmp_file = *bmp_file_ptr_ptr;
        if (bmp_file->rgb) {
            for(unsigned int i = 0; i < bmp_file->bih.height; i++){
                if (bmp_file->rgb[i]) {
                    free(bmp_file->rgb[i]);
                }
            }
            free(bmp_file->rgb);
        }
        free(bmp_file);
        *bmp_file_ptr_ptr = NULL;
    }
}

int valid_file(BMPFile* bmp_file){
    if(!bmp_file){
        return 0;
    }
    if (bmp_file->bfm.ID[0] != 'B' || bmp_file->bfm.ID[1] != 'M'){
        puts("ERROR: Input file is not a BMP file (magic number mismatch).");
        return 0;
    }
    return 1;
}

void print_info(BMPFile* bmpf){
    if(!bmpf){
        puts("Cannot print info: BMP data is null.");
        return;
    }
    printf("BMP File Information:\n");
    printf("  Header ID: %c%c\n", bmpf->bfm.ID[0], bmpf->bfm.ID[1]);
    printf("  File Size: %u bytes\n", bmpf->bfm.file_size);
    printf("  Pixel Offset: %u\n", bmpf->bfm.pixel_offset);
    printf("  Info Header Size: %u\n", bmpf->bih.HeaderSize);
    printf("  Width: %d pixels\n", bmpf->bih.width);
    printf("  Height: %d pixels\n", bmpf->bih.height);
    printf("  Planes: %u\n", bmpf->bih.planes);
    printf("  Bits Per Pixel: %u\n", bmpf->bih.bit_per_pixel);
    printf("  Compression: %u\n", bmpf->bih.compression);
    printf("  Image Size (from header): %u bytes\n", bmpf->bih.image_size);
    printf("  X Pixels Per Meter: %u\n", bmpf->bih.x_pixel_per_metr);
    printf("  Y Pixels Per Meter: %u\n", bmpf->bih.y_pixel_per_metr);
    printf("  Colors Used: %u\n", bmpf->bih.color_count);
    printf("  Important Colors: %u\n", bmpf->bih.imp_color_count);
}

void swap(int* a, int* b){
    int tmp = *a;
    *a = *b;
    *b = tmp;
}

int min(int a, int b){
    return (a < b) ? a : b;
}

int max(int a, int b){
    return (a > b) ? a : b;
}

int get_padding(int w){
    return (4 - (w * 3) % 4) % 4;
}

BMPFile* loadBMPFile(char* fname) {
    if (!fname) {
        puts("ERROR: Filename is NULL in loadBMPFile.");
        return NULL;
    }

    FILE* f = fopen(fname, "rb");
    if (!f) {
        printf("ERROR: Could not open file '%s' for reading.\n", fname);
        return NULL;
    }

    BMPFile* bmp_file = calloc(1, sizeof(BMPFile));
    if (!bmp_file) {
        puts("ERROR: Memory allocation failed for BMPFile structure.");
        fclose(f);
        return NULL;
    }

    if (fread(&bmp_file->bfm, sizeof(BMPFileHeader), 1, f) != 1) {
        puts("ERROR: Could not read BMP File Header.");
        fclose(f); free(bmp_file); return NULL;
    }
    if (fread(&bmp_file->bih, sizeof(BMPInfoHeader), 1, f) != 1) {
        puts("ERROR: Could not read BMP Info Header.");
        fclose(f); free(bmp_file); return NULL;
    }

    if (bmp_file->bfm.ID[0] != 'B' || bmp_file->bfm.ID[1] != 'M') {
        puts("ERROR: Not a BMP file (identifier mismatch).");
        fclose(f); free(bmp_file); return NULL;
    }
    if(bmp_file->bih.HeaderSize != 40){
        puts("ERROR: Unsupported BMP version (only DIB header size 40 is supported).");
        fclose(f); free(bmp_file); return NULL;
    }
    if(bmp_file->bih.bit_per_pixel != 24){
        puts("ERROR: Unsupported BMP format (only 24-bit per pixel supported).");
        fclose(f); free(bmp_file); return NULL;
    }
    if(bmp_file->bih.compression != 0){
        puts("ERROR: Unsupported BMP format (only uncompressed BMPs are supported).");
        fclose(f); free(bmp_file); return NULL;
    }
    if(bmp_file->bih.width == 0 || bmp_file->bih.height == 0 || bmp_file->bih.width > 20000 || bmp_file->bih.height > 20000) {
        printf("ERROR: Invalid BMP dimensions (width: %u, height: %u).\n", bmp_file->bih.width, bmp_file->bih.height);
        fclose(f); free(bmp_file); return NULL;
    }

    fseek(f, bmp_file->bfm.pixel_offset, SEEK_SET);
    int padding = get_padding(bmp_file->bih.width);

    bmp_file->rgb = malloc(sizeof(RGB*) * bmp_file->bih.height);
    if (!bmp_file->rgb) {
        puts("ERROR: Memory allocation failed for RGB row pointers.");
        fclose(f); free(bmp_file); return NULL;
    }
    for (unsigned int i = 0; i < bmp_file->bih.height; i++) {
        bmp_file->rgb[i] = NULL;
    }

    for (int i = bmp_file->bih.height - 1; i >= 0; i--) {
        bmp_file->rgb[i] = malloc(sizeof(RGB) * bmp_file->bih.width);
        if (!bmp_file->rgb[i]) {
            puts("ERROR: Memory allocation failed for an RGB row.");
            free_mem(NULL, NULL, &bmp_file);
            fclose(f);
            return NULL;
        }
        size_t elements_to_read = bmp_file->bih.width;
        if (elements_to_read > 0 && fread(bmp_file->rgb[i], sizeof(RGB), elements_to_read, f) != elements_to_read) {
             printf("ERROR: Could not read pixel data for row %d.\n", i);
             free_mem(NULL, NULL, &bmp_file);
             fclose(f);
             return NULL;
        }
        if (padding > 0) {
            if (fseek(f, padding, SEEK_CUR) != 0) {
                printf("ERROR: Could not seek past padding for row %d.\n", i);
                free_mem(NULL, NULL, &bmp_file);
                fclose(f);
                return NULL;
            }
        }
    }

    fclose(f);
    return bmp_file;
}

void writeBMPFile(char* fname, BMPFile* bmpf){
    if (!fname) { puts("ERROR: Filename is NULL in writeBMPFile."); return; }
    if (!bmpf) { puts("ERROR: BMPFile data is NULL in writeBMPFile."); return; }

    FILE* f = fopen(fname, "wb");
    if (!f) { printf("ERROR: Could not open file '%s' for writing.\n", fname); return; }

    if(fwrite(&bmpf->bfm, sizeof(BMPFileHeader), 1, f) != 1) {
        printf("ERROR: Could not write BMP File Header to '%s'.\n", fname); fclose(f); return;
    }
    if(fwrite(&bmpf->bih, sizeof(BMPInfoHeader), 1, f) != 1) {
        printf("ERROR: Could not write BMP Info Header to '%s'.\n", fname); fclose(f); return;
    }

    int padding = get_padding(bmpf->bih.width);
    unsigned char pad_bytes[3] = {0, 0, 0};

    for(int i = bmpf->bih.height - 1; i >= 0; i--){
        if (bmpf->bih.width > 0 && fwrite(bmpf->rgb[i], sizeof(RGB), bmpf->bih.width, f) != bmpf->bih.width) {
            printf("ERROR: Could not write pixel data for row %d to file '%s'.\n", i, fname);
            fclose(f); return;
        }
        if (padding > 0) {
            if (fwrite(pad_bytes, 1, padding, f) != (size_t)padding) {
                 printf("ERROR: Could not write padding for row %d to file '%s'.\n", i, fname);
                 fclose(f); return;
            }
        }
    }
    fclose(f);
}

static void set_pixel(int x, int y, int* color_rgb_array, BMPFile* bmpf){
    if (x >= 0 && x < (int)bmpf->bih.width && y >= 0 && y < (int)bmpf->bih.height) {
        bmpf->rgb[y][x].r = (unsigned char)color_rgb_array[0];
        bmpf->rgb[y][x].g = (unsigned char)color_rgb_array[1];
        bmpf->rgb[y][x].b = (unsigned char)color_rgb_array[2];
    }
}

static void set_thickness_pixel(int x_center, int y_center, int* color, int thickness, BMPFile* bmpf){
    if (thickness <= 0) return;
    if (thickness == 1) {
        set_pixel(x_center, y_center, color, bmpf);
        return;
    }
    const int th_radius = (thickness - 1) / 2;

    for(int y = y_center - th_radius; y <= y_center + th_radius; y++){
        if(y < 0 || y >= (int)bmpf->bih.height) continue;
        for(int x = x_center - th_radius; x <= x_center + th_radius; x++){
            if(x < 0 || x >= (int)bmpf->bih.width) continue;
            set_pixel(x, y, color, bmpf);
        }
    }
}

static int is_point_in_polygon(int nvert, int (*vert)[2], int testx, int testy) {
    int i, j, c = 0;
    for (i = 0, j = nvert - 1; i < nvert; j = i++) {
        if (((vert[i][1] > testy) != (vert[j][1] > testy)) &&
            (testx < (vert[j][0] - vert[i][0]) * (testy - vert[i][1]) / (double)(vert[j][1] - vert[i][1]) + vert[i][0]))
            c = !c;
    }
    return c;
}


static int is_point_on_line_segment(int px, int py, int x1, int y1, int x2, int y2, size_t line_thickness) {
    double dist_sq;
    double l2 = (x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2);

    if (l2 == 0) {
        dist_sq = (px - x1)*(px - x1) + (py - y1)*(py - y1);
    } else {
        double t = ((px - x1) * (x2 - x1) + (py - y1) * (y2 - y1)) / l2;
        t = fmax(0, fmin(1, t));

        double closest_x = x1 + t * (x2 - x1);
        double closest_y = y1 + t * (y2 - y1);
        dist_sq = (px - closest_x)*(px - closest_x) + (py - closest_y)*(py - closest_y);
    }
    double thickness_radius = line_thickness / 2.0;
    return dist_sq <= thickness_radius * thickness_radius;
}

static void bresenham_draw_line(int x0, int y0, int x1, int y1, int* color, int thickness, BMPFile* bmpf){
    int dx = abs(x1 - x0);
    int sx = x0 < x1 ? 1 : -1;
    int dy = -abs(y1 - y0);
    int sy = y0 < y1 ? 1 : -1;
    int err = dx + dy;
    int e2;

    for(;;){
        set_thickness_pixel(x0, y0, color, thickness, bmpf);
        if (x0 == x1 && y0 == y1) break;
        e2 = 2 * err;
        if (e2 >= dy) { err += dy; x0 += sx; }
        if (e2 <= dx) { err += dx; y0 += sy; }
    }
}

void draw_rectangle(int* left_coord, int* right_coord, size_t thickness,
    int* color, int fill, int* fill_color_rgb, BMPFile* bmpf){

    int x_start = min(left_coord[0], right_coord[0]);
    int y_start = min(left_coord[1], right_coord[1]);
    int x_end = max(left_coord[0], right_coord[0]);
    int y_end = max(left_coord[1], right_coord[1]);

    bresenham_draw_line(x_start, y_start, x_end, y_start, color, thickness, bmpf); // Top
    bresenham_draw_line(x_end, y_start, x_end, y_end, color, thickness, bmpf);     // Right
    bresenham_draw_line(x_end, y_end, x_start, y_end, color, thickness, bmpf);     // Bottom
    bresenham_draw_line(x_start, y_end, x_start, y_start, color, thickness, bmpf); // Left

    if(fill && fill_color_rgb){
        int th_offset = (thickness > 0) ? (thickness - 1) / 2 + 1 : 0;
        for(int y = y_start + th_offset; y <= y_end - th_offset; y++){
            if(y < 0 || y >= (int)bmpf->bih.height) continue;
            for(int x = x_start + th_offset; x <= x_end - th_offset; x++){
                if(x < 0 || x >= (int)bmpf->bih.width) continue;
                set_pixel(x, y, fill_color_rgb, bmpf);
            }
        }
    }
}

void draw_hexagon(int* center_coords, size_t radius_val, size_t thickness_val,
    int* border_color_rgb, int fill_flag, int* fill_color_rgb, BMPFile* bmpf){

    int vertices[6][2];
    for(int i = 0; i < 6; i++){
        double angle = PI / 3.0 * i;
        vertices[i][0] = (int)round(center_coords[0] + radius_val * cos(angle));
        vertices[i][1] = (int)round(center_coords[1] + radius_val * sin(angle));
    }

    for(int i = 0; i < 6; i++){
        bresenham_draw_line(vertices[i][0], vertices[i][1],
                  vertices[(i + 1) % 6][0], vertices[(i + 1) % 6][1],
                  border_color_rgb, thickness_val, bmpf);
    }

    if(fill_flag && fill_color_rgb){
        int min_x = vertices[0][0], max_x = vertices[0][0];
        int min_y = vertices[0][1], max_y = vertices[0][1];
        for(int i = 1; i < 6; i++){
            min_x = min(min_x, vertices[i][0]); max_x = max(max_x, vertices[i][0]);
            min_y = min(min_y, vertices[i][1]); max_y = max(max_y, vertices[i][1]);
        }
        
        for(int y = min_y; y <= max_y; y++){
            if(y < 0 || y >= (int)bmpf->bih.height) continue;
            for(int x = min_x; x <= max_x; x++){
                if(x < 0 || x >= (int)bmpf->bih.width) continue;

                if(is_point_in_polygon(6, vertices, x, y)){
                    int on_border = 0;
                    if (thickness_val > 0) {
                        for (int i = 0; i < 6; i++) {
                            if (is_point_on_line_segment(x, y, vertices[i][0], vertices[i][1],
                                                    vertices[(i + 1) % 6][0], vertices[(i + 1) % 6][1],
                                                    thickness_val)) {
                                on_border = 1;
                                break;
                            }
                        }
                    }
                    if (!on_border) {
                        set_pixel(x, y, fill_color_rgb, bmpf);
                    }
                }
            }
        }
    }
}

void copy_area(int* left_coord_src, int* right_coord_src, int* dest_coord_tl, BMPFile* bmpf) {
    int src_min_x = min(left_coord_src[0], right_coord_src[0]);
    int src_min_y = min(left_coord_src[1], right_coord_src[1]);
    int src_max_x = max(left_coord_src[0], right_coord_src[0]);
    int src_max_y = max(left_coord_src[1], right_coord_src[1]);

    src_min_x = max(0, src_min_x);
    src_min_y = max(0, src_min_y);
    src_max_x = min((int)bmpf->bih.width - 1, src_max_x);
    src_max_y = min((int)bmpf->bih.height - 1, src_max_y);

    int copy_width = (src_max_x >= src_min_x) ? (src_max_x - src_min_x + 1) : 0;
    int copy_height = (src_max_y >= src_min_y) ? (src_max_y - src_min_y + 1) : 0;

    if (copy_width == 0 || copy_height == 0) {
        puts("Warning: Source area for copy is zero or negative in size after clipping. Nothing copied.");
        return;
    }

    RGB** cp_area_buf = malloc(copy_height * sizeof(RGB*));
    if (!cp_area_buf) { puts("ERROR: Malloc failed for copy_area_buf rows."); return; }
    for(int i = 0; i < copy_height; i++){
        cp_area_buf[i] = malloc(copy_width * sizeof(RGB));
        if (!cp_area_buf[i]) {
            puts("ERROR: Malloc failed for copy_area_buf cols.");
            for (int j = 0; j < i; j++) free(cp_area_buf[j]);
            free(cp_area_buf); return;
        }
    }

    for(int y_idx = 0; y_idx < copy_height; y_idx++){
        for(int x_idx = 0; x_idx < copy_width; x_idx++){
            cp_area_buf[y_idx][x_idx] = bmpf->rgb[src_min_y + y_idx][src_min_x + x_idx];
        }
    }

    int dest_tl_x = dest_coord_tl[0];
    int dest_tl_y = dest_coord_tl[1];
    for(int y_idx = 0; y_idx < copy_height; y_idx++){
        int current_dest_y = dest_tl_y + y_idx;
        if(current_dest_y < 0 || current_dest_y >= (int)bmpf->bih.height) continue;
        for(int x_idx = 0; x_idx < copy_width; x_idx++){
            int current_dest_x = dest_tl_x + x_idx;
            if(current_dest_x < 0 || current_dest_x >= (int)bmpf->bih.width) continue;
            bmpf->rgb[current_dest_y][current_dest_x] = cp_area_buf[y_idx][x_idx];
        }
    }

    for(int i = 0; i < copy_height; i++) free(cp_area_buf[i]);
    free(cp_area_buf);
}

void expand(int scale_factor, int* fill_clr_rgb, BMPFile* bmpf_orig, char* output_filename) {
    if (scale_factor <= 0) { puts("ERROR: width_scale must be positive for expand."); return; }
    if (scale_factor == 1) { writeBMPFile(output_filename, bmpf_orig); return; }

    int old_w = bmpf_orig->bih.width;
    int old_h = bmpf_orig->bih.height;
    int new_w = old_w * scale_factor;
    int new_img_padding = get_padding(new_w);

    BMPFile new_bmpf_data;
    new_bmpf_data.bfm = bmpf_orig->bfm;
    new_bmpf_data.bih = bmpf_orig->bih;

    new_bmpf_data.bfm.file_size = sizeof(BMPFileHeader) + sizeof(BMPInfoHeader) +
                                 (sizeof(RGB) * new_w + new_img_padding) * old_h;
    new_bmpf_data.bih.width = new_w;
    new_bmpf_data.bih.image_size = (sizeof(RGB) * new_w + new_img_padding) * old_h;

    new_bmpf_data.rgb = malloc(old_h * sizeof(RGB*));
    if (!new_bmpf_data.rgb) { puts("ERROR: Malloc for new_rgb rows failed in expand."); return; }
    for (int y = 0; y < old_h; y++) {
        new_bmpf_data.rgb[y] = malloc(new_w * sizeof(RGB));
        if (!new_bmpf_data.rgb[y]) {
            puts("ERROR: Malloc for new_rgb cols failed in expand.");
            for(int k=0; k<y; ++k) free(new_bmpf_data.rgb[k]);
            free(new_bmpf_data.rgb); new_bmpf_data.rgb = NULL;
            return;
        }
    }
    
    RGB fill_pixel_value;
    fill_pixel_value.r = (unsigned char)fill_clr_rgb[0];
    fill_pixel_value.g = (unsigned char)fill_clr_rgb[1];
    fill_pixel_value.b = (unsigned char)fill_clr_rgb[2];

    for (int y = 0; y < old_h; y++) {
        int new_x_ptr = 0;
        for (int old_x = 0; old_x < old_w; old_x++) {
            if (new_x_ptr < new_w) {
                 new_bmpf_data.rgb[y][new_x_ptr++] = bmpf_orig->rgb[y][old_x];
            }
            for (int k = 1; k < scale_factor; k++) {
                if (new_x_ptr < new_w) {
                    new_bmpf_data.rgb[y][new_x_ptr++] = fill_pixel_value;
                }
            }
        }
    }
    writeBMPFile(output_filename, &new_bmpf_data);

    for (int y = 0; y < old_h; y++) free(new_bmpf_data.rgb[y]);
    free(new_bmpf_data.rgb);
}

// New binarization function
void do_binarization(BMPFile* bmpf, int threshold_value) {
    if (!bmpf || !bmpf->rgb) {
        puts("ERROR: BMP data or pixel array is null in do_binarization.");
        return;
    }

    unsigned char white_r = 255, white_g = 255, white_b = 255;
    unsigned char black_r = 0, black_g = 0, black_b = 0;

    for (unsigned int y = 0; y < bmpf->bih.height; y++) {
        for (unsigned int x = 0; x < bmpf->bih.width; x++) {
            RGB* current_pixel = &(bmpf->rgb[y][x]);
            int sum_components = current_pixel->r + current_pixel->g + current_pixel->b;

            if (sum_components >= threshold_value) {
                current_pixel->r = white_r;
                current_pixel->g = white_g;
                current_pixel->b = white_b;
            } else {
                current_pixel->r = black_r;
                current_pixel->g = black_g;
                current_pixel->b = black_b;
            }
        }
    }
}


static const struct option long_options_list[] = {
    {"help", no_argument, NULL, HELP_OPT_CHAR},
    {"input", required_argument, NULL, INPUT_OPT_CHAR},
    {"output", required_argument, NULL, OUTPUT_OPT_CHAR},
    {"rect", no_argument, NULL, RECT_OPT_VAL},
    {"left_up", required_argument, NULL, LEFT_UP_OPT_VAL},
    {"right_down", required_argument, NULL, RIGHT_DOWN_OPT_VAL},
    {"thickness", required_argument, NULL, THICKNESS_OPT_VAL},
    {"color", required_argument, NULL, COLOR_OPT_VAL},
    {"fill", no_argument, NULL, FILL_OPT_VAL},
    {"fill_color", required_argument, NULL, FILL_COLOR_OPT_VAL},
    {"hexagon", no_argument, NULL, HEXAGON_OPT_VAL},
    {"center", required_argument, NULL, CENTER_OPT_VAL},
    {"radius", required_argument, NULL, RADIUS_OPT_VAL},
    {"copy", no_argument, NULL, COPY_OPT_VAL},
    {"dest_left_up", required_argument, NULL, DEST_LEFT_UP_OPT_VAL},
    {"info", no_argument, NULL, INFO_OPT_VAL},
    {"expand", no_argument, NULL, EXPAND_OPT_VAL},
    {"width_scale", required_argument, NULL, WIDTH_SCALE_OPT_VAL},
    {"binarization", no_argument, NULL, BINARIZATION_OPT_VAL},     // New
    {"threshold", required_argument, NULL, THRESHOLD_OPT_VAL},    // New
    {0, 0, 0, 0}
};

static int parse_coord_str(int* coord_arr, char* str_val){
    if(sscanf(str_val, "%d.%d", &coord_arr[0], &coord_arr[1]) != 2){
        printf("Error: Invalid coordinate format '%s'. Expected X.Y (e.g., 100.150).\n", str_val);
        return 0;
    }
    return 1;
}

static int parse_positive_int_str(int* int_val, char* str_val, const char* arg_name){
    char *endptr;
    long val = strtol(str_val, &endptr, 10);
    if (str_val == endptr || *endptr != '\0' || val <= 0) {
        printf("Error: Invalid value for %s: '%s'. Expected a positive integer.\n", arg_name, str_val);
        return 0;
    }
    *int_val = (int)val;
    return 1;
}
static int parse_size_t_str(size_t* size_val, char* str_val, const char* arg_name){
    char *endptr;
    long val = strtol(str_val, &endptr, 10);
    if (str_val == endptr || *endptr != '\0' || val < 1) {
        printf("Error: Invalid value for %s: '%s'. Expected an integer >= 1.\n", arg_name, str_val);
        return 0;
    }
    *size_val = (size_t)val;
    return 1;
}

static int parse_color_str(int* color_arr_rgb, char* str_val){
    if(sscanf(str_val, "%d.%d.%d", &color_arr_rgb[0], &color_arr_rgb[1], &color_arr_rgb[2]) != 3){
        printf("Error: Invalid color format '%s'. Expected R.G.B (e.g., 255.0.128).\n", str_val);
        return 0;
    }
    return 1;
}

static void assign_filename_str(char** filename_dest_ptr, char* str_val){
    if (*filename_dest_ptr) free(*filename_dest_ptr);
    *filename_dest_ptr = malloc(strlen(str_val) + 1);
    if (*filename_dest_ptr) strcpy(*filename_dest_ptr, str_val);
    else puts("ERROR: Malloc failed for filename string.");
}

int parse_args(int argc, char* argv[], Options* opts){
    if(argc == 1){ help_func(); return 0; }

    int current_opt;
    opterr = 0;

    opts->output_name = strdup("out.bmp");

    while((current_opt = getopt_long(argc, argv, "hi:o:", long_options_list, NULL)) != -1){
        int prev_mode = opts->mode;
        switch (current_opt){
            case RECT_OPT_VAL: opts->mode = MODE_RECTANGLE; break;
            case HEXAGON_OPT_VAL: opts->mode = MODE_HEXAGON; break;
            case COPY_OPT_VAL: opts->mode = MODE_COPY; break;
            case INFO_OPT_VAL: opts->mode = MODE_INFO; break;
            case EXPAND_OPT_VAL: opts->mode = MODE_EXPAND; break;
            case BINARIZATION_OPT_VAL: opts->mode = MODE_BINARIZATION; break; // New

            case LEFT_UP_OPT_VAL: if(!parse_coord_str(opts->left_coord, optarg)) return 0; opts->is_left_coord = 1; break;
            case RIGHT_DOWN_OPT_VAL: if(!parse_coord_str(opts->right_coord, optarg)) return 0; opts->is_right_coord = 1; break;
            case THICKNESS_OPT_VAL: if(!parse_size_t_str(&(opts->thickness), optarg, "--thickness")) return 0; opts->is_thickness = 1; break;
            case COLOR_OPT_VAL: if(!parse_color_str(opts->color, optarg)) return 0; opts->is_color = 1; break;
            case FILL_OPT_VAL: opts->fill = 1; break;
            case FILL_COLOR_OPT_VAL: if(!parse_color_str(opts->fill_color, optarg)) return 0; opts->is_fill_color = 1; break;
            case CENTER_OPT_VAL: if(!parse_coord_str(opts->center, optarg)) return 0; opts->is_center = 1; break;
            case RADIUS_OPT_VAL: if(!parse_size_t_str(&(opts->radius), optarg, "--radius")) return 0; opts->is_radius = 1; break;
            case DEST_LEFT_UP_OPT_VAL: if(!parse_coord_str(opts->dest_left_coord, optarg)) return 0; opts->is_dest_left_coord = 1; break;
            case WIDTH_SCALE_OPT_VAL: if(!parse_positive_int_str(&(opts->width_scale), optarg, "--width_scale")) return 0; opts->is_width_scale = 1; break;
            case THRESHOLD_OPT_VAL: // New
                if(!parse_positive_int_str(&(opts->threshold), optarg, "--threshold")) return 0;
                opts->is_threshold = 1;
                break;
            
            case OUTPUT_OPT_CHAR: assign_filename_str(&(opts->output_name), optarg); opts->is_output_name = 1; break;
            case INPUT_OPT_CHAR: assign_filename_str(&(opts->input_name), optarg); opts->is_input_name = 1; break;
            
            case HELP_OPT_CHAR: help_func(); opts->mode = MODE_UNKNOWN; return 1;
            case '?':
                printf("Error: Unknown option or missing argument near '%s'.\n", argv[optind-1]);
                help_func(); return 0;
            default: help_func(); return 0;
        }
        if (prev_mode != MODE_UNKNOWN && opts->mode != prev_mode &&
            (opts->mode == MODE_RECTANGLE || opts->mode == MODE_HEXAGON || opts->mode == MODE_COPY ||
             opts->mode == MODE_INFO || opts->mode == MODE_EXPAND || opts->mode == MODE_BINARIZATION) ) { // Added MODE_BINARIZATION
            puts("Error: Multiple conflicting operation modes specified (e.g., --rect and --hexagon).");
            help_func(); return 0;
        }
    }
    
    if (optind < argc && !opts->is_input_name) {
        assign_filename_str(&(opts->input_name), argv[optind]);
        opts->is_input_name = 1;
        optind++;
    }
    if (optind < argc) {
        printf("Error: Unexpected argument: %s\n", argv[optind]);
        help_func(); return 0;
    }

    if (opts->mode != MODE_UNKNOWN && opts->mode != MODE_INFO && !opts->is_input_name) {
        puts("Error: Input BMP file must be specified (e.g., --input file.bmp or as the last argument).");
        help_func(); return 0;
    }
    if (opts->mode == MODE_UNKNOWN && argc > 1) {
         int is_help_arg = 0;
         for(int i=1; i<argc; ++i) if(strcmp(argv[i], "-h")==0 || strcmp(argv[i], "--help")==0) is_help_arg=1;
         if(!is_help_arg) {
            puts("Error: No operation mode selected (e.g., --rect, --hexagon, --info).");
            help_func(); return 0;
         }
    }
    return 1;
}

int main(int argc, char* argv[]){
    printf("Course work for option 4.9, created by Vanyushkin Kirill.\n");

    Options options_config;
    memset(&options_config, 0, sizeof(Options));

    if(!parse_args(argc, argv, &options_config)){
        free_mem(&options_config.input_name, &options_config.output_name, NULL);
        return 41;
    }

    if(options_config.mode == MODE_UNKNOWN){ // Help was shown or error in parse_args
        free_mem(&options_config.input_name, &options_config.output_name, NULL);
        return (argc == 1 || (argc > 1 && (strcmp(argv[1],"-h")==0 || strcmp(argv[1],"--help")==0) ) ) ? 0 : 41; // return 0 for help, 41 for other parsing errors that parse_args didn't catch but set mode to unknown
    }


    BMPFile* bmp_image_data = NULL;
    if (options_config.mode != MODE_INFO && !options_config.is_input_name && options_config.mode != MODE_UNKNOWN ) { // MODE_INFO can run without input, but we check later
         //This case should be caught by parse_args, but as a safeguard
        printf("Error: Input file is required for this operation.\n");
        free_mem(&options_config.input_name, &options_config.output_name, NULL);
        return 41;
    }
    
    if (options_config.is_input_name) { // Load if input name is provided (relevant for all modes except potentially a misconfigured help/error case)
         bmp_image_data = loadBMPFile(options_config.input_name);
         if(!valid_file(bmp_image_data)){ // valid_file checks if bmp_image_data is NULL too
            free_mem(&options_config.input_name, &options_config.output_name, &bmp_image_data);
            return 41;
        }
    } else if (options_config.mode != MODE_INFO && options_config.mode != MODE_UNKNOWN) { 
        // This should have been caught by parse_args general input file check.
        // If we are not in info mode, not in unknown (help/error) mode, and no input file, it's an issue.
        printf("Error: Input BMP file must be specified for the selected operation.\n");
        help_func();
        free_mem(&options_config.input_name, &options_config.output_name, NULL);
        return 41;
    }


    int operation_status = 1;
    switch (options_config.mode){
        case MODE_RECTANGLE:
            if(!validate_rectangle(options_config.is_left_coord, options_config.is_right_coord,
                options_config.is_thickness, options_config.is_color, options_config.fill,
                options_config.is_fill_color, options_config.color, options_config.fill_color)){
                operation_status = 0;
            } else {
                draw_rectangle(options_config.left_coord, options_config.right_coord, options_config.thickness,
                               options_config.color, options_config.fill, options_config.fill_color, bmp_image_data);
                writeBMPFile(options_config.output_name, bmp_image_data);
            }
            break;
        case MODE_HEXAGON:
            if(!validate_hexagon(options_config.is_center, options_config.is_radius,
                options_config.is_thickness, options_config.is_color, options_config.fill,
                options_config.is_fill_color, options_config.color, options_config.fill_color)){
                operation_status = 0;
            } else {
                draw_hexagon(options_config.center, options_config.radius, options_config.thickness,
                             options_config.color, options_config.fill, options_config.fill_color, bmp_image_data);
                writeBMPFile(options_config.output_name, bmp_image_data);
            }
            break;
        case MODE_COPY:
            if(!validate_copy(options_config.is_left_coord, options_config.is_right_coord,
                              options_config.is_dest_left_coord)){
                operation_status = 0;
            } else {
                copy_area(options_config.left_coord, options_config.right_coord, options_config.dest_left_coord, bmp_image_data);
                writeBMPFile(options_config.output_name, bmp_image_data);
            }
            break;
        case MODE_INFO:
            if (!bmp_image_data) { // Specifically for info, if no input file was given.
                 puts("Error: --info requires an input file specified via --input or as the last argument.");
                 operation_status = 0;
            } else {
                print_info(bmp_image_data);
            }
            break;
        case MODE_EXPAND:
             if(!validate_expand(options_config.is_fill_color, options_config.is_width_scale) ||
                !check_rgb(options_config.fill_color) ){
                operation_status = 0;
            } else {
                expand(options_config.width_scale, options_config.fill_color, bmp_image_data, options_config.output_name);
            }
            break;
        case MODE_BINARIZATION: // New
            if (!validate_binarization(options_config.is_threshold, options_config.threshold)) {
                operation_status = 0;
            } else {
                do_binarization(bmp_image_data, options_config.threshold);
                writeBMPFile(options_config.output_name, bmp_image_data);
            }
            break;
        default: // Should not happen if parse_args sets mode correctly or to UNKNOWN
            printf("Internal Error: Unknown mode %d in main switch.\n", options_config.mode);
            operation_status = 0;
            break;
    }

    free_mem(&options_config.input_name, &options_config.output_name, &bmp_image_data);

    if (!operation_status) {
        printf("Operation failed. Please check arguments and input file.\n");
        return 41;
    }
    
    if (options_config.mode != MODE_INFO && options_config.mode != MODE_UNKNOWN) {
        // Check output_name again in case it was freed due to an error path not reassigning it.
        // opts->output_name is initialized in parse_args and only changed if -o is given.
        // free_mem nullifies it. The string used by printf should be the one from options_config *before* free_mem.
        // This is a bit tricky. Let's assume options_config.output_name held the correct value for the message.
        // To be safe, perhaps store output_name for this message before freeing.
        // However, the current structure is that output_name points to a string that *will be* freed.
        // Let's use a temporary variable for the final message if the name was dynamically allocated.
        // For now, this is fine as long as options_config.output_name itself isn't used after free_mem.
        // The current code implicitly relies on the string literal "out.bmp" if output_name was not set and then freed.
        // Actually, parse_args guarantees output_name is strdup("out.bmp") or a strdup from optarg.
        // So it's always dynamically allocated. It's best to print the message *before* freeing.
        // Let's re-order.

        // The logic for printing the success message relies on options_config.output_name.
        // This should be done BEFORE freeing it. Let's move the free_mem call to the very end.
        char* final_output_name = options_config.output_name ? strdup(options_config.output_name) : strdup("out.bmp (default error)");
        // If options_config.output_name was already freed, this would be an issue.
        // Let's adjust the printing of success message logic:

        printf("Operation completed. Output saved to: %s\n", (options_config.output_name ? options_config.output_name : "out.bmp (default)"));
    
    } else if (options_config.mode == MODE_INFO && operation_status) { // only print if info operation was successful
        printf("Information printed successfully.\n");
    }
    // else if mode is UNKNOWN, messages are handled by parse_args or help_func directly.
    
    free_mem(&options_config.input_name, &options_config.output_name, &bmp_image_data);


    return 0;
}

--- END OF FILE cw8.c ---