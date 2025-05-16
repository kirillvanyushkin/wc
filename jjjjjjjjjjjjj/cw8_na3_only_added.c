// В структуре Options:
    int threshold; // For binarization
    int is_threshold; // Flag for binarization threshold

// Прототипы новых функций:
int validate_binarization(int is_threshold, int threshold_val); // New declaration
void do_binarization(BMPFile* bmpf, int threshold_value); // New declaration

// Новые #define константы:
#define BINARIZATION_OPT_VAL 1016 // New
#define THRESHOLD_OPT_VAL 1017    // New
#define MODE_BINARIZATION 6 // New

// В функции help_func():
    puts("  --binarization              Binarize the image (set pixels to black or white)."); // New
    puts("    --threshold N             Threshold value for binarization (integer > 0). Sum of R+G+B >= N -> white, else black."); // New
    puts("  ./cw --input image.bmp --binarization --threshold 384 -o binarized_image.bmp"); // New example

// Новая функция validate_binarization():
// New validation function
int validate_binarization(int is_threshold, int threshold_val){
    int success = 1;
    if(!is_threshold) { puts("Error: --binarization requires --threshold."); success = 0; }
    // The check for threshold_val > 0 is already handled by parse_positive_int_str
    return success;
}

// Новая функция do_binarization():
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

// В массиве long_options_list[]:
    {"binarization", no_argument, NULL, BINARIZATION_OPT_VAL},     // New
    {"threshold", required_argument, NULL, THRESHOLD_OPT_VAL},    // New

// В функции parse_args(), внутри switch(current_opt):
            case BINARIZATION_OPT_VAL: opts->mode = MODE_BINARIZATION; break; // New

            case THRESHOLD_OPT_VAL: // New
                if(!parse_positive_int_str(&(opts->threshold), optarg, "--threshold")) return 0;
                opts->is_threshold = 1;
                break;

// В функции parse_args(), в проверке конфликта режимов:
             opts->mode == MODE_INFO || opts->mode == MODE_EXPAND || opts->mode == MODE_BINARIZATION) ) { // Added MODE_BINARIZATION

// В функции main(), внутри switch(options_config.mode):
        case MODE_BINARIZATION: // New
            if (!validate_binarization(options_config.is_threshold, options_config.threshold)) {
                operation_status = 0;
            } else {
                do_binarization(bmp_image_data, options_config.threshold);
                writeBMPFile(options_config.output_name, bmp_image_data);
            }
            break;