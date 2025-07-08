#include <stdio.h>
#include <stdlib.h>

/** Geospatial Data Abstraction Library (GDAL) **/
#include "gdal.h"       // public (C callable) GDAL entry points
#include "cpl_conv.h"   // various convenience functions for CPL
#include "cpl_string.h" // various convenience functions for strings


#include "utils/const.h"
#include "utils/alloc.h"
#include "utils/dir.h"
#include "utils/string.h"



void usage(char *exe, int exit_code){

  printf("\n");
  printf("Usage: %s *files\n", exe);
  printf("  \n");
  printf("  *files can be one or multiple input files of the same dimensions\n");
  printf("  The last band is used for maximum-X compositing\n");
  printf("  Commonly, it is maximum-NDVI\n");

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
  int n_input;
  char **input_path;
  char output_path[STRLEN];
} args_t;


void parse_args(int argc, char *argv[], args_t *args){
int opt;

  opterr = 0;

  while ((opt = getopt(argc, argv, "o:")) != -1){
    switch(opt){
      case 'o':
        copy_string(args->output_path, STRLEN, optarg);
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


  args->n_input = argc-optind;

  if (args->n_input == 0) {
    fprintf(stderr, "no input files specified\n");
    usage(argv[0], FAILURE);
  }

  alloc_2D((void***)&args->input_path, args->n_input, STRLEN, sizeof(char));

  for (int i = 0; i < args->n_input; i++) {
    copy_string(args->input_path[i], STRLEN, argv[optind + i]);
    if (!fileexist(args->input_path[i])) {
      fprintf(stderr, "file %s does not exist\n", args->input_path[i]);
      usage(argv[0], FAILURE);
    }
  }

  return;
}


int main ( int argc, char *argv[] ){


args_t args;

  parse_args(argc, argv, &args);

  GDALAllRegister();

  image_t *images = NULL;
  alloc((void**)&images, args.n_input, sizeof(image_t));


  for (int i = 0; i < args.n_input; i++) {

    GDALDatasetH dataset;

    if ((dataset = GDALOpen(args.input_path[i], GA_ReadOnly)) == NULL){ 
      fprintf(stderr, "could not open %s\n", args.input_path[i]); 
      usage(argv[0], FAILURE);
    }

    images[i].ncol  = GDALGetRasterXSize(dataset);
    images[i].nrow  = GDALGetRasterYSize(dataset);
    images[i].ncell = images[i].ncol*images[i].nrow;
    
    copy_string(images[i].projection, STRLEN, GDALGetProjectionRef(dataset));
    GDALGetGeoTransform(dataset, images[i].geotransformation);


    images[i].nband = GDALGetRasterCount(dataset);
    alloc_2D((void***)&images[i].image, images[i].nband, images[i].ncell, sizeof(short));

    for (int b = 0; b < images[i].nband; b++) {

      GDALRasterBandH band;

      band = GDALGetRasterBand(dataset, b+1);
      //int has_nodata = 0;

      //images[i].nodata = GDALGetRasterNoDataValue(band, &has_nodata);
      //if (!has_nodata) {
      //  fprintf(stderr, "%s has no nodata value.\n", args.input_path[i]); 
      //  usage(argv[0], FAILURE);
      //}
  
      images[i].datatype = GDALGetRasterDataType(band);
      if (images[i].datatype > GDT_Int16) {
        printf("datatype needs to ne Int16 (is: %s)\n", GDALGetDataTypeName(images[i].datatype)); 
        usage(argv[0], FAILURE);
      }

      if (GDALRasterIO(band, GF_Read, 0, 0, images[i].ncol, images[i].nrow, images[i].image[b], 
          images[i].ncol, images[i].nrow, GDT_Int16, 0, 0) == CE_Failure){
        printf("could not read band %d from %s\n", b+1, args.input_path[i]); 
        usage(argv[0], FAILURE);
      }

    }

    printf("file: %s\n", args.input_path[i]);
    printf("projection: %s\n", images[i].projection);
    printf("origin: %.6f %.6f\n", images[i].geotransformation[0], images[i].geotransformation[3]);
    printf("resolution: %.6f %.6f\n", images[i].geotransformation[1], images[i].geotransformation[5]);
    printf("dimensions: %d x %d = %d pixels\n", images[i].nrow, images[i].ncol, images[i].ncell);
    printf("bands: %d\n", images[i].nband);
    //printf("nodata: %f\n", images[i].nodata);
    printf("datatype: %s\n", GDALGetDataTypeName(images[i].datatype));
    printf("\n");

    GDALClose(dataset);

  }


  if (args.n_input > 1) {

    for (int i = 1; i < args.n_input; i++) {

      if (images[i].ncol != images[0].ncol || 
          images[i].nrow != images[0].nrow ||
          images[i].nband != images[0].nband) {
        fprintf(stderr, "input files have different dimensions\n");
        usage(argv[0], FAILURE);
      }

    }
  
  }


  short **composite = NULL;
  alloc_2D((void***)&composite, images[0].nband, images[0].ncell, sizeof(short));

  for (int c = 0; c < images[0].ncell; c++) {

    short maximum = SHRT_MIN;
    int i_maximum = -1;
    
    for (int i = 0; i < args.n_input; i++) {
      
      int skip = 0;
      //printf("processing cell %d of %d (%s)\n", c+1, images[0].ncell, args.input_path[i]);
      //printf("  maximum: %d\n", maximum);
      //printf("  i_maximum: %d\n", i_maximum); 
      //printf("  current: %d\n", images[i].image[images[i].nband-1][c]);
      for (int b = 0; b < images[0].nband-1; b++) {
        if (images[i].image[b][c] == SHRT_MIN ||
            images[i].image[b][c] == SHRT_MAX) {
          skip = 1;
          break;
        }
      }

      if (skip) continue;

      if (images[i].image[images[i].nband-1][c] != 0 &&
          images[i].image[images[i].nband-1][c] > maximum) {
          maximum = images[i].image[images[i].nband-1][c];
          i_maximum = i;
      }
        
    }

    //printf("found best\n");
    //printf("  maximum: %d\n", maximum);
    //printf("  i_maximum: %d\n", i_maximum); 

    for (int b = 0; b < images[0].nband; b++) {

      //printf("  band %d: ", b+1);
      //printf("  reflectance: %d\n", images[i_maximum].image[b][c]);

      if (i_maximum < 0) {
        composite[b][c] = SHRT_MIN;
      } else {
        composite[b][c] = images[i_maximum].image[b][c];
      }

    }

    if (composite[0][c] == 32767) printf("issue in cell %d\n", c); 

  }




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


  if ((output_dataset = GDALCreate(output_driver, args.output_path, images[0].ncol, images[0].nrow, images[0].nband-1, images[0].datatype, output_options)) == NULL) {
    printf("Error creating file %s.\n", args.output_path);
    usage(argv[0], FAILURE);
  }

  for (int b = 0; b < images[0].nband-1; b++) {

    output_band = GDALGetRasterBand(output_dataset, b+1);
    GDALSetRasterNoDataValue(output_band, SHRT_MIN);

    if (GDALRasterIO(output_band, GF_Write, 0, 0, images[0].ncol, images[0].nrow, 
      composite[b], images[0].ncol, images[0].nrow, GDT_Int16, 0, 0) == CE_Failure){
      printf("Unable to write band %d in %s.\n", b, args.output_path); 
      usage(argv[0], FAILURE);
    }

  }


  GDALSetGeoTransform(output_dataset, images[0].geotransformation);
  GDALSetProjection(output_dataset,   images[0].projection);

  GDALClose(output_dataset);


  for (int i = 0; i < args.n_input; i++) {
    free_2D((void**)images[i].image, images[i].nband);
  }
  free((void*)images);

  if (output_options != NULL) CSLDestroy(output_options);   

  return SUCCESS;

}
