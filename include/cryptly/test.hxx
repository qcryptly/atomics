#pragma once

#include <gtest/gtest.h>

class Test : public ::testing::Test {};

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
