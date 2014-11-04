#pragma once
#include <vector>
#include <string>

namespace
{
  // finds the least significant bit, returns ~0 if number is 0
  unsigned int find_first_set(ulong val) {
    // remove all bits above lowest bit
    val &= -val;

    unsigned int index = -1;
    while (val != 0)
    {
      val >>= 1;
      index++;
    }
    return index;
  }
}

class gf2
{
  unsigned int rows, cols, ulong_width;
  std::vector<std::vector<ulong>> matrix;
  public:
    gf2(unsigned int rows, unsigned int cols)
      : rows(rows), cols(cols), ulong_width(sizeof(ulong) * 8),
        matrix(cols, std::vector<ulong>(rows, 0))
    {
    }

    bool get_bit(unsigned int row, unsigned int col)
    {
      return (matrix[col][row / ulong_width] & ((ulong)1 << row % ulong_width)) > 0;
    }

    void add_bit(unsigned int row, unsigned int col)
    {
      matrix[col][row / ulong_width] ^= (ulong)1 << row % ulong_width;
    }

    // col1 = col1 + col2
    void add_columns(unsigned int col1, unsigned int col2)
    {
      for (unsigned int r = 0; r <= rows / ulong_width; r++)
      {
        matrix[col1][r] ^= matrix[col2][r];
      }
    }

    std::vector<unsigned int> fast_gauss()
    {
      std::vector<unsigned int> marks;
      for (unsigned int j = 0; j < cols; j++)
      {
        for (unsigned int row = 0; row <= rows / ulong_width; row++)
        {
          if (matrix[j][row] == 0)
          {
            continue;
          }
          unsigned int i = find_first_set(matrix[j][row]);
          marks.push_back(i);

          for (unsigned int k = 0; k < cols; k++)
          {
            if (k == j)
            {
              continue;
            }

            if (get_bit(i, k) == true)
            {
              add_columns(k, j);
            }
          }
          break;
        }
      }

      return marks;
    }

    std::string to_string()
    {
      std::string res;

      for (unsigned int row = 0; row < rows; row++)
      {
        for (unsigned int col = 0; col < cols; col++)
        {
          res += get_bit(row, col) ? "1 " : "0 ";
        }
        res += "\n";
      }

      return res;
    }
};
