#include <algorithm>
#include <cassert>
#include <iostream>
#include <map>
#include <memory>
#include <numeric>
#include <queue>
#include <set>
#include <stack>
#include <type_traits>
#include <utility>
#include <vector>

template <typename... Ts>
struct is_listkind_helper
{
};

template <typename T, typename _ = void>
struct is_listkind : public std::false_type
{
};

template <typename T>
struct is_listkind<T,
                   std::conditional_t<
                       false,
                       is_listkind_helper<
                           typename T::value_type,
                           typename T::size_type,
                           typename T::allocator_type,
                           typename T::iterator,
                           typename T::const_iterator,
                           decltype(std::declval<T>().size()),
                           decltype(std::declval<T>().begin()),
                           decltype(std::declval<T>().end()),
                           decltype(std::declval<T>().cbegin()),
                           decltype(std::declval<T>().cend())>,
                       void>> : public std::true_type
{
};

template <typename T, typename _ = void>
struct Source
{
    std::set<T> alphabets;
    std::map<T, double> freqs;

    Source<std::vector<T>> extend(size_t n)
    {
        Source<std::vector<T>> res;
        if (n == 0)
            return res;

        size_t depth = 0;
        std::stack<std::pair<std::vector<T>, double>> stack;

        for (auto &&e : alphabets)
        {
            stack.emplace({e}, freqs.at(e));
        }

        std::set<std::vector<T>> reached;
        while (!stack.empty())
        {
            auto t = stack.top();
            auto w = t.first;
            auto p = t.second;
            if (w.size() < n)
            {
                for (auto &&e : alphabets)
                {
                    w.push_back(e);
                    stack.emplace(w, p * freqs.at(e));
                    w.pop_back();
                }
            }
            if (w.size() == n)
            {
                res.alphabets.emplace(t);
            }
            stack.pop();
        }
    }
};

template <typename T>
struct Source<is_listkind_helper<
    typename T::value_type,
    typename T::size_type,
    typename T::allocator_type,
    typename T::iterator,
    typename T::const_iterator,
    decltype(std::declval<T>().size()),
    decltype(std::declval<T>().begin()),
    decltype(std::declval<T>().end()),
    decltype(std::declval<T>().cbegin()),
    decltype(std::declval<T>().cend())>>
{
    std::set<T> alphabets;
    std::map<T, double> freqs;

    Source<T> extend(size_t n)
    {
        Source<T> res;
        if (n == 0)
            return res;

        size_t depth = 0;
        std::stack<std::pair<std::vector<T>, double>> stack;

        std::set<std::vector<T>> reached;
        while (!stack.empty())
        {
            auto t = stack.top();
            auto w = t.first;
            auto p = t.second;
            if (w.size() < n)
            {
                for (auto &&e : alphabets)
                {
                    w.push_back(e);
                    stack.emplace(w, p * freqs.at(e));
                    w.pop_back();
                }
            }
            if (w.size() == n)
            {
                res.alphabets.emplace(t);
            }
            stack.pop();
        }
    }
};

size_t weight(const std::vector<bool> &code)
{
    return std::accumulate(code.begin(), code.end(), 0, [](size_t acc, bool e) -> size_t
                           { return acc + static_cast<size_t>(e); });
}

size_t dist_hamming(const std::vector<bool> &c1, const std::vector<bool> &c2)
{
    size_t res = 0;
    for (size_t i = 0; i < c1.size() && i < c2.size(); i++)
    {
        res += static_cast<size_t>(c1.at(i) and c2.at(i));
    }
    return res;
}

template <typename T>
class Huffman
{
private:
    struct Compare
    {
        bool operator()(const std::pair<std::set<T>, size_t> &a, const std::pair<std::set<T>, size_t> &b) const
        {
            return a.second > b.second;
        }
    };

public:
    std::map<T, std::vector<bool>> code;
    Huffman(const Source<T> &s)
    {
        std::priority_queue<std::pair<std::set<T>, double>, std::vector<std::pair<std::set<T>, double>>, Compare> pq;
        for (auto &&e : s.alphabets)
        {
            pq.push(std::make_pair(std::set{e}, s.freqs.at(e)));
            code.emplace(e, std::vector<bool>{});
        }
        while (pq.size() > 1)
        {
            auto pq1 = pq.top();
            pq.pop();
            auto pq2 = pq.top();
            pq.pop();
            for (auto &&e : pq1.first)
            {
                code.at(e).insert(code.at(e).begin(), 0);
            }
            for (auto &&e : pq2.first)
            {
                code.at(e).insert(code.at(e).begin(), 1);
            }
            std::set<T> uni = pq1.first;
            uni.insert(pq2.first.begin(), pq2.first.end());
            pq.push(std::make_pair(uni, pq1.second + pq2.second));
        }
    }
};

template <typename T, T _operator(const T &, const T &), T _inverse(const T &), T _identity>
struct Group
{
    using type = T;

    static const T identity = _identity;

    static T append(const T &x, const T &y)
    {
        return _operator(x, y);
    }

    static T inverse(const T &x)
    {
        return _inverse(x);
    }

    static T exp(const T &x, size_t n)
    {
        T res = x;
        for (size_t i = 1; i < n; i++)
        {
            res = append(res, x);
        }
        return res;
    }

    static T concat(const std::vector<T> &xs)
    {
        T res = xs.front();
        for (size_t i = 1; i < xs.size(); i++)
        {
            res = append(res, xs.at(i));
        }
        return res;
    }
};

template <typename T, T _add(const T &, const T &), T _neg(const T &), T _add_id, T _mul(const T &, const T &), T _rec(const T &), T _mul_id>
struct Field
{
    using type = T;
    static T add(const T &a, const T &b)
    {
        return _add(a, b);
    }
    static T sub(const T &a, const T &b)
    {
        return _add(a, _neg(b));
    }
    static T mul(const T &a, const T &b)
    {
        return _mul(a, b);
    }
    static T div(const T &a, const T &b)
    {
        return _mul(a, _rec(b));
    }
    static T neg(const T &a)
    {
        return _neg(a);
    }
    ///reciprocal
    static T rec(const T &a)
    {
        return _rec(a);
    }
    static Group<T, add, _neg, _add_id> additive;
    static Group<T, mul, _rec, _mul_id> multiplicative;
    static const T add_id = _add_id;
    static const T mul_id = _mul_id;
};

template <typename K>
class Matrix
{
private:
    using T = K::type;
    std::vector<T> data;
    K field;

public:
    const size_t line, column;
    Matrix(size_t n) : line(n), column(n)
    {
        for (size_t i = 0; i < n; i++)
        {
            for (size_t j = 0; j < n; j++)
            {
                data.push_back((i == j ? K::multiplicative.identity : K::additive.identity));
            }
        }
    }
    Matrix(size_t line, size_t column, const std::initializer_list<std::initializer_list<T>> &data = {}) : line(line), column(column)
    {
        for (auto &&line : data)
        {
            for (auto &&e : line)
            {
                this->data.push_back(e);
            }
        }
        if (data.size() == 0)
        {
            this->data.resize(line * column);
        }
    }

    Matrix(const std::vector<T> &v) : data(v), line(1), column(v.size())
    {
    }

    operator std::vector<T>()
    {
        return data;
    }

    std::vector<T>::const_reference at(size_t i, size_t j) const
    {
        return data.at(i * column + j);
    }

    std::vector<T>::reference at(size_t i, size_t j)
    {
        return data.at(i * column + j);
    }

    Matrix<K> scalared(const T &k)
    {
        Matrix<K> res(line, column);
        for (size_t i = 0; i < line; i++)
        {
            for (size_t j = 0; j < column; j++)
            {
                res.at(i, j) = K::mul(at(i, j), k);
            }
        }
        return res;
    }

    Matrix<K> operator*(const Matrix<K> &m) const
    {
        Matrix<K> res(line, m.column);
        for (size_t i = 0; i < line; i++)
        {
            for (size_t j = 0; j < m.column; j++)
            {
                for (size_t k = 0; k < column; k++)
                {
                    res.at(i, j) = K::add(res.at(i, j), K::mul(at(i, k), m.at(k, j)));
                }
            }
        }
        return res;
    }

    Matrix<K> transposed() const
    {
        Matrix<K> res(column, line);
        for (size_t i = 0; i < column; i++)
        {
            for (size_t j = 0; j < line; j++)
            {
                res.at(i, j) = at(j, i);
            }
        }
        return res;
    }

    void swap_lines(size_t i1, size_t i2)
    {
        for (size_t j = 0; j < column; j++)
        {
            std::swap(at(i1, j), at(i2, j));
        }
    }

    void scale_line(size_t i, const T &l)
    {
        for (size_t j = 0; j < column; j++)
        {
            at(i, j) = K::mul(at(i, j), l);
        }
    }

    void add_scaledline2line(size_t i1, const T &k, size_t i2)
    {
        for (size_t j = 0; j < column; j++)
        {
            at(i2, j) = K::add(K::mul(at(i1, j), k), at(i2, j));
        }
    }

    void gj_eliminate()
    {
        size_t k = 0;
        for (size_t j = 0; j < column; j++)
        {
            size_t i = k;
            for (; i < line && at(i, j) == K::additive.identity; i++)
            {
            }
            if (i < line)
            {
                if (i != k)
                    swap_lines(k, i);
                scale_line(k, K::multiplicative.inverse(at(k, j)));
                for (i = 0; i < line; i++)
                {
                    if (i != k)
                        add_scaledline2line(k, K::multiplicative.inverse(at(i, j)), i);
                }
                k++;
            }
        }
    }

    Matrix<K> at_line(size_t i)
    {
        Matrix<K> l(1, column);
        for (size_t j = 0; j < column; j++)
        {
            l.at(i, j) = at(i, j);
        }
        return l;
    }

    Matrix<K> at_column(size_t j)
    {
        Matrix<K> c(line, 1);
        for (size_t i = 0; i < line; i++)
        {
            c.at(i, j) = at(i, j);
        }
        return c;
    }

    Matrix<K> concat_lines(const Matrix<K> &mat)
    {
        Matrix<K> res(line + mat.line, column);
        for (size_t i = 0; i < res.line; i++)
        {
            for (size_t j = 0; j < res.column; j++)
            {
                if (i < line)
                    res.at(i, j) = at(i, j);
                else
                    res.at(i, j) = mat.at(i - line, j);
            }
        }
        return res;
    }

    Matrix<K> concat_columns(const Matrix<K> &mat)
    {
        Matrix<K> res(line, column + mat.column);
        for (size_t i = 0; i < res.line; i++)
        {
            for (size_t j = 0; j < res.column; j++)
            {
                if (j < column)
                    res.at(i, j) = at(i, j);
                else
                    res.at(i, j) = mat.at(i, j - column);
            }
        }
        return res;
    }

    Matrix<K> select_lines(size_t begin, size_t end)
    {
        Matrix<K> res(end - begin + 1, column);
        for (size_t i = begin; i <= end; i++)
        {
            for (size_t j = 0; j < column; j++)
            {
                res.at(i - begin, j) = at(i, j);
            }
        }
        return res;
    }

    Matrix<K> select_columns(size_t begin, size_t end)
    {
        Matrix<K> res(line, end - begin + 1);
        for (size_t i = 0; i < line; i++)
        {
            for (size_t j = begin; j <= end; j++)
            {
                res.at(i, j - begin) = at(i, j);
            }
        }
        return res;
    }

    std::pair<Matrix<K>, Matrix<K>> cut_square()
    {
        if (line < column)
        {
            return std::make_pair(select_columns(0, line - 1), select_columns(line, column - 1));
        }
        return std::make_pair(select_lines(0, column - 1), select_lines(column, line));
    }

    Matrix<K> gj_eliminated() const
    {
        Matrix<K> res = *this;
        res.gj_eliminate();
        return res;
    }

    void print() const
    {
        for (size_t i = 0; i < line; i++)
        {
            for (size_t j = 0; j < column; j++)
            {
                std::cout << at(i, j) << " ";
            }
            std::cout << "\n";
        }
    }
};

template <typename K>
class Code
{
private:
public:
    using type = K::type;
    const Matrix<K> gen, parity;
    const size_t n, k;
    Code(const Matrix<K> &gen) : gen(gen.gj_eliminated()),
                                 parity(gen.gj_eliminated().cut_square().second.transposed().concat_columns(Matrix<K>(gen.column - gen.line))),
                                 n(gen.column), k(gen.line)
    {
    }

    std::vector<type> code(const std::vector<type> &x)
    {
        return Matrix<K>(x) * gen;
    }

    std::vector<type> check(const std::vector<type> &v)
    {
        return Matrix<K>(v) * parity.transposed();
    }
};

bool logic_add(const bool &a, const bool &b)
{
    return a xor b;
}

bool logic_mul(const bool &a, const bool &b)
{
    return a and b;
}

bool logic_neg(const bool &a)
{
    return not a;
}

bool logic_id(const bool &a)
{
    return a;
}

using GF2 = Field<bool, logic_add, logic_neg, false, logic_mul, logic_id, true>;

int main()
{
    Source<char> s;
    s.alphabets.insert('A');
    s.alphabets.insert('B');
    s.alphabets.insert('C');
    s.alphabets.insert('D');
    s.freqs['A'] = 5;
    s.freqs['B'] = 4;
    s.freqs['C'] = 2;
    s.freqs['D'] = 2;
    Huffman hf(s);

    Matrix<GF2> mat(3, 3, {
                              {1, 1, 0},
                              {1, 0, 1},
                              {1, 1, 1},
                          });
    mat.gj_eliminate();

    Matrix<GF2> G(4, 7, {
                            {1, 0, 0, 0, 0, 1, 1},
                            {0, 1, 0, 0, 1, 0, 1},
                            {0, 0, 1, 0, 1, 1, 0},
                            {0, 0, 0, 1, 1, 1, 1},
                        });
    Code hamming74(G);
    hamming74.parity.print();
    std::vector<bool> x = {0, 1, 0, 1};
    auto v = hamming74.code(x);
    std::cout << "\n";
    Matrix<GF2>(v).print();
    std::cout << "\n";
    Matrix<GF2>(hamming74.check(v)).print();
    return 0;
}