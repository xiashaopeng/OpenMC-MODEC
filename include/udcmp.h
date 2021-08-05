#ifndef UDCMP_H_
#define UDCMP_H_

#include <iostream>
#include <functional>
#include <vector>
#include <iterator>

template<typename T>
class udcmp : public std::binary_function<T, T, bool> {
    std::vector<bool>& v_;
    bool b_;
public:
    explicit udcmp(std::vector<bool>& v) : v_(v) {}
    bool operator () (const T& t1, const T& t2) {
        b_ = t1 < t2;
        v_.push_back(b_);
        return b_;
    }
};

// ref"https://stackoverflow.com/questions/37368787/c-sort-one-vector-based-on-another-one"
namespace sort_helper {

    template <typename _Data, typename _Order>
    struct value_reference_t;

    template <typename _Data, typename _Order>
    struct value_t {
        _Data data;
        _Order val;
        inline value_t(_Data _data, _Order _val) : data(_data), val(_val) {}
        inline value_t(const value_reference_t<_Data, _Order>& rhs);
    };

    template <typename _Data, typename _Order>
    struct value_reference_t {
        typename std::vector<_Data>::iterator pdata;
        typename std::vector<_Order>::iterator pval;
        value_reference_t(typename std::vector<_Data>::iterator _itData, typename std::vector<_Order>::iterator _itVal) : pdata(_itData), pval(_itVal) {}
        inline value_reference_t& operator = (const value_reference_t& rhs) { *pdata = *rhs.pdata; *pval = *rhs.pval; return *this; }
        inline value_reference_t& operator = (const value_t<_Data, _Order>& rhs) { *pdata = rhs.data; *pval = rhs.val; return *this; }
        inline bool operator < (const value_reference_t& rhs) { return *pval < *rhs.pval; }
    };

    template <typename _Data, typename _Order>
    struct value_iterator_t :
        std::iterator< std::random_access_iterator_tag, value_t<_Data, _Order>, ptrdiff_t, value_t<_Data, _Order>*, value_reference_t<_Data, _Order> >
    {
        typename std::vector<_Data>::iterator itData;
        typename std::vector<_Order>::iterator itVal;
        value_iterator_t(typename std::vector<_Data>::iterator _itData, typename std::vector<_Order>::iterator _itVal) : itData(_itData), itVal(_itVal) {}
        inline ptrdiff_t operator - (const value_iterator_t& rhs) const { return itVal - rhs.itVal; }
        inline value_iterator_t operator + (ptrdiff_t off) const { return value_iterator_t(itData + off, itVal + off); }
        inline value_iterator_t operator - (ptrdiff_t off) const { return value_iterator_t(itData - off, itVal - off); }
        inline value_iterator_t& operator ++ () { ++itData; ++itVal; return *this; }
        inline value_iterator_t& operator -- () { --itData; --itVal; return *this; }
        inline value_iterator_t operator ++ (int) { return value_iterator_t(itData++, itVal++); }
        inline value_iterator_t operator -- (int) { return value_iterator_t(itData--, itVal--); }
        inline value_t<_Data, _Order> operator * () const { return value_t<_Data, _Order>(*itData, *itVal); }
        inline value_reference_t<_Data, _Order> operator * () { return value_reference_t<_Data, _Order>(itData, itVal); }
        inline bool operator  < (const value_iterator_t& rhs) const { return itVal < rhs.itVal; }
        inline bool operator == (const value_iterator_t& rhs) const { return itVal == rhs.itVal; }
        inline bool operator != (const value_iterator_t& rhs) const { return itVal != rhs.itVal; }
    };

    template <typename _Data, typename _Order>
    inline value_t<_Data, _Order>::value_t(const value_reference_t<_Data, _Order>& rhs)
        : data(*rhs.pdata), val(*rhs.pval) {}

    template <typename _Data, typename _Order>
    bool operator < (const value_t<_Data, _Order>& lhs, const value_reference_t<_Data, _Order>& rhs) {
        return lhs.val < *rhs.pval;
    }

    template <typename _Data, typename _Order>
    bool operator < (const value_reference_t<_Data, _Order>& lhs, const value_t<_Data, _Order>& rhs) {
        return *lhs.pval < rhs.val;
    }

    template <typename _Data, typename _Order>
    void swap(value_reference_t<_Data, _Order> lhs, value_reference_t<_Data, _Order> rhs) {
        std::swap(*lhs.pdata, *rhs.pdata);
        std::swap(*lhs.pval, *rhs.pval);
    }

} // namespace sort_helper


#endif