#include <algorithm>
#include <array>
#include <iostream>

template <int Cols>
struct RowVec {
    std::array<float, Cols> row;

    float& operator [] (int n) { return row[n]; }
    float operator [] (int n) const { return row[n]; }

    float x() { static_assert(Cols > 0, "out of bounds"); return row[0]; }
    float y() { static_assert(Cols > 1, "out of bounds"); return row[1]; }
    float z() { static_assert(Cols > 2, "out of bounds"); return row[2]; }

    RowVec& operator /= (float rhs)
    {
        for (int n = 0; n < Cols; n++)
            row[n] /= rhs;
        return *this;
    }

    RowVec& operator -= (const RowVec& rhs)
    {
        for (int n = 0; n < Cols; n++)
            row[n] -= rhs[n];
        return *this;
    }
};

template <int Cols>
RowVec<Cols> operator * (float lhs, RowVec<Cols> rhs)
{
    for (int n = 0; n < Cols; n++)
        rhs[n] *= lhs;
    return rhs;
}

template <int Cols>
std::ostream& operator << (std::ostream& oss, const RowVec<Cols> row)
{
    oss << "[" << row[0];
    for (int n = 1; n < Cols; n++)
        oss << ", " << row[n];
    return oss << "]";
}

template <int Rows>
struct ColVec {
    std::array<float, Rows> col;

    float& operator [] (int m) { return col[m]; }
    float operator [] (int m) const { return col[m]; }

    float x() { static_assert(Rows > 0, "index out of bounds"); return col[0]; }
    float y() { static_assert(Rows > 1, "index out of bounds"); return col[1]; }
    float z() { static_assert(Rows > 2, "index out of bounds"); return col[2]; }
};

template <int Rows>
ColVec<Rows> operator + (ColVec<Rows> lhs, const ColVec<Rows>& rhs)
{
    for (int m = 0; m < Rows; m++)
        lhs[m] += rhs[m];
    return lhs;
}

template <int Rows>
ColVec<Rows> operator * (float lhs, ColVec<Rows> rhs)
{
    for (int m = 0; m < Rows; m++)
        rhs[m] *= lhs;
    return rhs;
}

template <int Rows>
std::ostream& operator << (std::ostream& oss, const ColVec<Rows> col)
{
    oss << "<" << col[0];
    for (int m = 1; m < Rows; m++)
        oss << ", " << col[m];
    return oss << ">";
}

template <int Rows, int Cols>
struct Mat {
    std::array<RowVec<Cols>, Rows> rows;

    RowVec<Cols>& operator [] (int m) { return rows[m]; }
    const RowVec<Cols>& operator [] (int m) const { return rows[m]; }

    const ColVec<Rows> col(int n) const
    {
        ColVec<Rows> column;
        for (int m = 0; m < Rows; m++)
            column[m] = rows[m][n];
        return column;
    }
};

template <int Rows, int Cols>
std::ostream& operator << (std::ostream& oss, const Mat<Rows, Cols> mat)
{
    oss << "[" << mat[0];
    for (int m = 1; m < Rows; m++)
        oss << ", " << mat[m];
    return oss << "]";
}

template <int N>
float operator * (const RowVec<N>& lhs, const ColVec<N>& rhs)
{
    float result = 0;
    for (int i = 0; i < N; i++)
        result += lhs[i] * rhs[i];
    return result;
}

template <int Rows, int Shared>
ColVec<Rows> operator * (const Mat<Rows, Shared>& lhs, const ColVec<Shared>& rhs)
{
    ColVec<Rows> result;
    for (int m = 0; m < Rows; m++)
        result[m] = lhs[m] * rhs;
    return result;
}

template <int Shared, int Cols>
RowVec<Cols> operator * (const RowVec<Shared>& lhs, const Mat<Shared, Cols>& rhs)
{
    RowVec<Cols> result;
    for (int n = 0; n < Cols; n++)
        result[n] = lhs * rhs.col(n);
    return result;
}

template <int Rows, int Shared, int Cols>
Mat<Rows, Cols> operator * (const Mat<Rows, Shared>& lhs, const Mat<Shared, Cols>& rhs)
{
    Mat<Rows, Cols> result;
    for (int m = 0; m < Rows; m++)
        for (int n = 0; n < Cols; n++)
            result[m][n] = lhs[m] * rhs.col(n);
    return result;
}

Mat<3, 3> translate(float x, float y)
{
    return {
        1.f, 0.f,   x,
        0.f, 1.f,   y,
        0.f, 0.f, 1.f,
    };
}

Mat<3, 3> scale(float x, float y)
{
    return {
          x, 0.f, 0.f,
        0.f,   y, 0.f,
        0.f, 0.f, 1.f,
    };
}

/**
Computes the inverse of a square matrix using Gaussian eliminations.
*/
template <int N>
Mat<N, N> inverse(const Mat<N, N>& input)
{
    // Augment input with an identity matrix on the right
    Mat<N, N * 2> block;
    for (int m = 0; m < N; m++)
        for (int n = 0; n < N * 2; n++)
        {
            if (n < N)
                block[m][n] = input[m][n];
            else
                block[m][n] = m == n - N;
        }

    // Implementation of
    // https://en.wikipedia.org/wiki/Row_echelon_form#Pseudocode_for_reduced_row_echelon_form
    int lead = 0;
    const int row_count = N;
    const int column_count = N * 2;
    for (int r = 0; r < row_count; r++)
    {
        if (column_count <= lead)
            break;

        int i = r;
        while (block[i][lead] == 0)
        {
            i = i + 1;
            if (row_count == i)
            {
                i = r;
                lead = lead + 1;
                if (column_count == lead)
                    break;
            }
        }
        if (i != r)
            std::swap(block[i], block[r]);

        block[r] /= block[r][lead];
        for (int j = 0; j < row_count; j++)
            if (j != r)
                block[j] -= block[j][lead] * block[r];
        lead = lead + 1;
    }

    Mat<N, N> output;
    for (int m = 0; m < N; m++)
        for (int n = 0; n < N; n++)
            output[m][n] = block[m][n + N];
    return output;
}
