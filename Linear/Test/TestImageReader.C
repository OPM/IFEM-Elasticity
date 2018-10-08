//==============================================================================
//!
//! \file TestImageReader.C
//!
//! \date Sept 2018
//!
//! \author Kjetil A. Johannessen / SINTEF
//!
//! \brief Tests for image reader
//!
//==============================================================================

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#include "gtest/gtest.h"


TEST(TestImageReader, png)
{
  int width, height, nrChannels;
  unsigned char *image = stbi_load("refdata/chess.png",
                                   &width,
                                   &height,
                                   &nrChannels, 0);
  ASSERT_EQ(width,      512);
  ASSERT_EQ(height,     512);
  ASSERT_EQ(nrChannels,   4);
  int k = 0;
  for(int j=0; j<height; j++) {
    for(int i=0; i<width; i++) {
      bool is_white = ((i/64)%2 == (j/64)%2);
      ASSERT_EQ(image[k  ], is_white * 255); // red
      ASSERT_EQ(image[k+1], is_white * 255); // green
      ASSERT_EQ(image[k+2], is_white * 255); // blue
      ASSERT_EQ(image[k+3], 255);            // alpha
      k += 4;
    }
  }
}

TEST(TestImageReader, rectangular)
{
  int width, height, nrChannels;
  unsigned char *image = stbi_load("refdata/rectangular_chess.png",
                                   &width,
                                   &height,
                                   &nrChannels, 0);
  ASSERT_EQ(width,      256);
  ASSERT_EQ(height,     128);
  ASSERT_EQ(nrChannels,   4);
  int k = 0;
  for(int j=0; j<height; j++) {
    for(int i=0; i<width; i++) {
      bool is_white = ((i/64)%2 == (j/64)%2);
      ASSERT_EQ(image[k  ], is_white * 255); // red
      ASSERT_EQ(image[k+1], is_white * 255); // green
      ASSERT_EQ(image[k+2], is_white * 255); // blue
      ASSERT_EQ(image[k+3], 255);            // alpha
      k += 4;
    }
  }
}
