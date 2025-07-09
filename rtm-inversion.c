#include <stdio.h>
#include <stdlib.h>
#include <float.h>


/** Geospatial Data Abstraction Library (GDAL) **/
#include "gdal.h"       // public (C callable) GDAL entry points
#include "cpl_conv.h"   // various convenience functions for CPL
#include "cpl_string.h" // various convenience functions for strings


#include "utils/const.h"
#include "utils/alloc.h"
#include "utils/dir.h"
#include "utils/string.h"
#include "utils/table.h"

//#include <omp.h>

void usage(char *exe, int exit_code){

  printf("\n");
  printf("Usage: %s -l LUT.csv -s simulations.csv -i input.tif -o output.tif [-a 0.01] [-n 100]\n", exe);
  printf("  \n");
  printf("  adapt file names\n");
  printf("  -a inversion stops when accuracy is met\n");
  printf("  -n inversion stops when max iterations are used\n");
  printf("   use -a 0 to disable accuracy check, this brute-forces the inversion\n");
  printf("\n");

  exit(exit_code);
  return;
}

typedef struct {
  GDALDataType datatype;
  int nrow, ncol, ncell, nband;
  char projection[STRLEN];
  double geotransformation[6];
  //double nodata;
  short **image;
} image_t;

typedef struct {
  char lut_path[STRLEN];
  char simulation_path[STRLEN];
  char input_path[STRLEN];
  char output_path[STRLEN];
  int max_iterations;
  float accuracy;
} args_t;


void parse_args(int argc, char *argv[], args_t *args){
int opt;

  opterr = 0;

  copy_string(args->lut_path, STRLEN, "NULL");
  copy_string(args->simulation_path, STRLEN, "NULL");
  copy_string(args->input_path, STRLEN, "NULL");
  copy_string(args->output_path, STRLEN, "NULL");

  args->accuracy = 0.01;
  args->max_iterations = 100;

  while ((opt = getopt(argc, argv, "l:s:i:o:a:n:")) != -1){
    switch(opt){
      case 'l':
        copy_string(args->lut_path, STRLEN, optarg);
        break;
      case 's':
        copy_string(args->simulation_path, STRLEN, optarg);
        break;
      case 'i':
        copy_string(args->input_path, STRLEN, optarg);
        break;
      case 'o':
        copy_string(args->output_path, STRLEN, optarg);
        break;
      case 'a':
        args->accuracy = atof(optarg);
        break;
      case 'n':
        args->max_iterations = atoi(optarg);
        break;
      case '?':
        if (isprint(optopt)){
          fprintf(stderr, "Unknown option `-%c'.\n", optopt);
        } else {
          fprintf(stderr, "Unknown option character `\\x%x'.\n", optopt);
        }
        usage(argv[0], FAILURE);
      default:
        fprintf(stderr, "Error parsing arguments.\n");
        usage(argv[0], FAILURE);
    }
  }

  if (strcmp(args->lut_path, "NULL") == 0 ||
      strcmp(args->simulation_path, "NULL") == 0 ||
      strcmp(args->input_path, "NULL") == 0 ||
      strcmp(args->output_path, "NULL") == 0) {
    fprintf(stderr, "missing arguments\n");
    usage(argv[0], FAILURE);
  }

  int n_input = argc-optind;

  if (n_input > 0) {
    fprintf(stderr, "too many parameters specified\n");
    usage(argv[0], FAILURE);
  }

  return;
}


int main ( int argc, char *argv[] ){


args_t args;

  parse_args(argc, argv, &args);

  table_t lut = read_table(args.lut_path, false, false);
  table_t simulations = read_table(args.simulation_path, false, false);

  if (lut.nrow != simulations.nrow) {
    fprintf(stderr, "LUT and simulations have different number of rows (%d vs %d)\n", lut.nrow, simulations.nrow);
    usage(argv[0], FAILURE);
  }

  for (int i = 0; i < simulations.nrow; i++) {
    for (int j = 0; j < simulations.ncol; j++) {
      simulations.data[i][j] *= 10000;
    }
  }


  print_table(&lut, true, false);
  print_table(&simulations, true, false);

  //for (int i = 0; i < simulations.nrow; i++) {
  //  for (int j = 0; j < simulations.ncol; j++) {
  //    simulations.data[i][j] = (simulations.data[i][j] - simulations.mean[j]) / simulations.sd[j];
  //  }
  //}

  GDALAllRegister();

  GDALDatasetH dataset;
  image_t input;

  if ((dataset = GDALOpen(args.input_path, GA_ReadOnly)) == NULL){ 
    fprintf(stderr, "could not open %s\n", args.input_path); 
    usage(argv[0], FAILURE);
  }

  input.ncol  = GDALGetRasterXSize(dataset);
  input.nrow  = GDALGetRasterYSize(dataset);
  input.ncell = input.ncol*input.nrow;
  
  copy_string(input.projection, STRLEN, GDALGetProjectionRef(dataset));
  GDALGetGeoTransform(dataset, input.geotransformation);


  input.nband = GDALGetRasterCount(dataset);

  if (input.nband != simulations.ncol) {
    fprintf(stderr, "number of bands (%d) does not match number of simulations (%d)\n", input.nband, simulations.ncol);
    usage(argv[0], FAILURE);
  }

  alloc_2D((void***)&input.image, input.nband, input.ncell, sizeof(short));

  for (int b = 0; b < input.nband; b++) {

    GDALRasterBandH band;

    band = GDALGetRasterBand(dataset, b+1);
    //int has_nodata = 0;

    //input.nodata = GDALGetRasterNoDataValue(band, &has_nodata);
    //if (!has_nodata) {
    //  fprintf(stderr, "%s has no nodata value.\n", args.input_path); 
    //  usage(argv[0], FAILURE);
    //}

    input.datatype = GDALGetRasterDataType(band);
    if (input.datatype > GDT_Int16) {
      printf("datatype needs to ne Int16 (is: %s)\n", GDALGetDataTypeName(input.datatype)); 
      usage(argv[0], FAILURE);
    }

    if (GDALRasterIO(band, GF_Read, 0, 0, input.ncol, input.nrow, input.image[b], 
        input.ncol, input.nrow, GDT_Int16, 0, 0) == CE_Failure){
      printf("could not read band %d from %s\n", b+1, args.input_path); 
      usage(argv[0], FAILURE);
    }

  }

  printf("file: %s\n", args.input_path);
  printf("projection: %s\n", input.projection);
  printf("origin: %.6f %.6f\n", input.geotransformation[0], input.geotransformation[3]);
  printf("resolution: %.6f %.6f\n", input.geotransformation[1], input.geotransformation[5]);
  printf("dimensions: %d x %d = %d pixels\n", input.nrow, input.ncol, input.ncell);
  printf("bands: %d\n", input.nband);
  //printf("nodata: %f\n", input.nodata);
  printf("datatype: %s\n", GDALGetDataTypeName(input.datatype));
  printf("\n");

  GDALClose(dataset);





  float **inversion = NULL;
  alloc_2D((void***)&inversion, lut.ncol+1, input.ncell, sizeof(float));

  
  srand(time(NULL)); // initialize random seed
  
  //#pragma omp parallel shared(input, inversion, lut, simulations, args) default(none)
  //{
  //unsigned int seed = time(NULL) ^ omp_get_thread_num();
// Use rand_r(&seed) for random numbers in this thread
//  #pragma omp for
  for (int c = 0; c < input.ncell; c++) {

    for (int o = 0; o < lut.ncol; o++) {
      inversion[o][c] = -1.0;
    }
    inversion[lut.ncol][c] = -1.0; // store -1.0 for mae in addtional band


    int skip = 0;

    for (int b = 0; b < input.nband; b++) {
      if (input.image[b][c] == SHRT_MIN ||
          input.image[b][c] == SHRT_MAX) {
        skip = 1;
        break;
      }
    }

    if (skip) continue;
    

    float min_mae = FLT_MAX;
    int i_min_mae = -1;
    int ctr = 0;

    // brute-force inversion
    if (args.accuracy <= FLT_EPSILON) {

      for (int i = 0; i < simulations.nrow; i++) {

        float mae = 0.0;

        for (int b = 0; b < input.nband; b++) {
          //float ref = (input.image[b][c] - simulations.mean[b]) / simulations.sd[b];
          float ref = (float)input.image[b][c];
          mae += fabs(ref - simulations.data[i][b]);
        }

        mae /= input.nband;

        if (mae < min_mae) {
          min_mae = mae;
          i_min_mae = i;
        }

      }

    // use accuracy to early-stop inversion
    } else {

      while (min_mae > args.accuracy && ctr < args.max_iterations) {

        // randomly select a row from the LUT
        int i = rand() % simulations.nrow;
        
        float mae = 0.0;

        for (int b = 0; b < input.nband; b++) {

          //float ref = (input.image[b][c] - simulations.mean[b]) / simulations.sd[b];
          float ref = (float)input.image[b][c];
          mae += fabs(ref - simulations.data[i][b]);
          
        }

        mae /= input.nband;

        if (mae < min_mae) {
          min_mae = mae;
          i_min_mae = i;
        }

        ctr++;

      }

    }

    //printf("cell %d: min mae = %.2f at row %d. %d iterations used\n", c, min_mae, i_min_mae, ctr);

    if (i_min_mae >= 0) {
      for (int o = 0; o < lut.ncol; o++) {
        inversion[o][c] = lut.data[i_min_mae][o];
      }
      inversion[lut.ncol][c] = min_mae; // store min mae in addtional band
    }

  }

  //}

  GDALDatasetH output_dataset = NULL;
  GDALRasterBandH output_band = NULL;
  GDALDriverH output_driver = NULL;
  char **output_options = NULL;

  if ((output_driver = GDALGetDriverByName("GTiff")) == NULL) {
    printf("%s driver not found\n", "GTiff"); 
    usage(argv[0], FAILURE);
  }

  output_options = CSLSetNameValue(output_options, "COMPRESS", "ZSTD");
  output_options = CSLSetNameValue(output_options, "PREDICTOR", "2");
  output_options = CSLSetNameValue(output_options, "BIGTIFF", "YES");
  //output_options = CSLSetNameValue(output_options, "OVERVIEWS", "NONE");


  if ((output_dataset = GDALCreate(output_driver, args.output_path, input.ncol, input.nrow, lut.ncol+1, GDT_Float32, output_options)) == NULL) {
    printf("Error creating file %s.\n", args.output_path);
    usage(argv[0], FAILURE);
  }

  for (int o = 0; o < lut.ncol+1; o++) {

    output_band = GDALGetRasterBand(output_dataset, o+1);
    GDALSetRasterNoDataValue(output_band, -1.0);

    if (GDALRasterIO(output_band, GF_Write, 0, 0, input.ncol, input.nrow, 
      inversion[o], input.ncol, input.nrow, GDT_Float32, 0, 0) == CE_Failure){
      printf("Unable to write band %d in %s.\n", o+1, args.output_path); 
      usage(argv[0], FAILURE);
    }

  }

  // one additional band for the mae
  output_band = GDALGetRasterBand(output_dataset, lut.ncol+1);
  GDALSetRasterNoDataValue(output_band, -1.0);

  if (GDALRasterIO(output_band, GF_Write, 0, 0, input.ncol, input.nrow, 
    inversion[lut.ncol], input.ncol, input.nrow, GDT_Float32, 0, 0) == CE_Failure){
    printf("Unable to write band %d in %s.\n", lut.ncol+1, args.output_path); 
    usage(argv[0], FAILURE);
  }


  GDALSetGeoTransform(output_dataset, input.geotransformation);
  GDALSetProjection(output_dataset,   input.projection);

  GDALClose(output_dataset);

  free_2D((void**)input.image, input.nband);
  free_2D((void**)inversion, lut.ncol+1);

  free_table(&lut);
  free_table(&simulations);
  
  if (output_options != NULL) CSLDestroy(output_options);   

  return SUCCESS;

}
