typedef unsigned long ulong;

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

// based on "A Fast Algorithm for Gaussian Elimination over GF(2)"
// https://www.cs.umd.edu/~gasarch/TOPICS/factoring/fastgauss.pdf

class gf2
{
  unsigned int rows, cols, ulong_width;
  std::vector<std::vector<ulong>> matrix;
  public:
    gf2(unsigned int rows, unsigned int cols)
      : rows(rows), cols(cols), ulong_width(sizeof(ulong) * 8),
        matrix(cols, std::vector<ulong>((rows - 1) / ulong_width + 1, 0))
    {
    }

    bool get_bit(unsigned int row, unsigned int col)
    {
      return (matrix[col][row / ulong_width] & (1L << row % ulong_width)) > 0;
    }

    void add_bit(unsigned int row, unsigned int col)
    {
      matrix[col][row / ulong_width] ^= 1L << row % ulong_width;
    }

    // col1 = col1 + col2
    void add_columns(unsigned int col1, unsigned int col2)
    {
      for (unsigned int r = 0; r <= rows / ulong_width; r++)
      {
        matrix[col1][r] ^= matrix[col2][r];
      }
    }

    std::vector<std::vector<unsigned int>> fast_gauss()
    {
      auto dependencies = std::vector<unsigned int>(cols, ~0);
      auto marks = std::set<unsigned int>();

      for (unsigned int j = 0; j < cols; j++)
      {
        for (unsigned int row = 0; row <= rows / ulong_width; row++)
        {
          if (matrix[j][row] == 0)
          {
            continue;
          }
          unsigned int i = find_first_set(matrix[j][row]) + row * ulong_width;
          marks.insert(i);
          // keep track of dependencies so we don't have to recalculate them later
          dependencies[j] = i;

          for (unsigned int k = 0; k < cols; k++)
          {
            if (k == j)
            {
              continue;
            }

            if (get_bit(i, k))
            {
              add_columns(k, j);
            }
          }
          break;
        }
      }
      auto result = std::vector<std::vector<unsigned int>>();

      for (unsigned int row = 0; row < rows; row++)
      {
        if (marks.count(row) != 0)
        {
          continue;
        }

        auto res = std::vector<unsigned int>();
        res.push_back(row);
        for (unsigned int col = 0; col < cols; col++)
        {
          if (get_bit(row, col))
          {
            res.push_back(dependencies[col]);
          }
        }

        result.push_back(res);
      }

      return result;
    }

    std::string to_string()
    {
      std::string res;
      res += std::to_string(rows) + " " + std::to_string(cols) + "\n";

      for (unsigned int row = 0; row < rows; row++)
      {
        for (unsigned int col = 0; col < cols; col++)
        {
          res += get_bit(row, col) ? "1" : "0";
        }
        res += "\n";
      }

      return res;
    }
};
