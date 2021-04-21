#define MTL_WITH_RANGEDFOR 1
#define MTL_WITH_INITLIST 1
#define MTL_WITH_AUTO 1
#define MTL_WITH_MOVE 1
#define MTL_WITH_STATICASSERT 1
#define MTL_WITH_TEMPLATE_ALIAS 1
#define MTL_WITH_VARIADIC_TEMPLATE 1
//#define MTL_WITH_DEFAULTIMPL 1

#include <cmath>
#include <iostream>
#include <list>
#include <iterator>
#include <vector>
#include <algorithm> // find
#include <cmath>
#include <tuple>
#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/itl/itl.hpp>
#include <boost/numeric/mtl/utility/tag.hpp>

using namespace mtl;
using std::cout;
using std::endl;
using std::ostream;
using std::prev;
using std::next;
using std::distance;
using std::copy;
using std::ostream_iterator;
using std::vector;
using std::find;
using std::begin;
using std::end;
using std::pow;

void calc_l(int n, int m)
{
    int l= 0;
    for (int i= 1; i < n; i++){
        for (int j= 1; j < m; j++) {
            l= i + (m - 1 - j)*(n-1);
            cout << "l(" << i << "," << j << ") = " << l << endl;
        }
    }
}

void calc_lr(int n, int m)
{
    int l= 0;
    for ( int j= m-1; j >= 1; j-- ){
        for ( int i= 1; i < n; i++ ) { 
            l= i + (m - 1 - j)*(n-1) - 1;
            cout << "l(" << i << "," << j << ") = " << l << endl;
        }
    }
}
/*
def calc_l(n,m):
...     for i in range(1,n):
...         for j in range(1,m):
...             l = i + (m-1-j)*(n-1)
...             print('l({:d},{:d}) = {:d}'.format(i,j,l))
*/

template <typename T>
class Point2D
{
    public:
    T x, y;
    Point2D(T a, T b) : x(a), y(b) {}
};

template <typename T>
class PointND
{
    public:
    int length;
    mtl::dense_vector<T> value;
    PointND(std::initializer_list<T> alist) : length(alist.size())
    {
        //std::copy(alist.begin(), alist.end(), value);
        value = alist;
        if (length == 4) {
            std::cout << "PointND size is 4. We are in Space-Time Continuum!\n";
        }
        else {
            cout << "size is: " << (value) << endl;
        }
    }
};

/*
template <typename T>
std::ostream& operator<<(std::ostream& os, const std::list< std::list<T> >& alist)
// Return os that can use type Point2D, and pass to std::cout
{
    for ( auto it= alist.begin(), e= alist.end(); it < e; ++it)
        os << "{ " << *it[0] << " , " << *it[1] << " }" << endl;
    return os;
}
*/

template <typename T>
std::ostream& operator<<(std::ostream& os, const Point2D<T>& pt)
/* Return os that can use type Point2D, and pass to std::cout  */
{
    os << "(" << pt.x << "," << pt.y << ")" << endl;
    return os;
}

template <typename T>
std::ostream& operator<<(std::ostream& os, const PointND<T>& pt)
/* Return os that can use type PointND, and pass to std::cout  */
{
    os << "(";
    //for ( auto xi : pt.value ) {
    for ( auto it= std::begin(pt.value), e= std::end(pt.value); it < e; ++it ) {
        if ( it < e - 1)
            os << *it << ",";
        else
            os << *it << ")\n";
    }

    return os;
}

//---------------------------------------------//
// Print Vectors<T> and Vector<Vector<int>>
//---------------------------------------------//
template <typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& v)
{
    os << "[";
    for ( auto vi= std::begin(v), ve= std::end(v); vi < ve; ++vi ) {
        if ( vi < ve - 1)
            os << *vi << " , ";
        else
            os << *vi << "]\n";
    }
    return os;
}

template <typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<std::vector<T> >& vv)
{
    // print vector or vector of vectors
    for ( auto it= std::begin(vv), e= std::end(vv); it < e; ++it ) {
        os << "[";
        for ( auto vi= std::begin(*it), ve= std::end(*it); vi < ve; ++vi ) {
            if ( vi < ve - 1)
                os <<  *vi << " , ";
            else
                os << *vi << "]\n";
        }
    }
return os;
}

template <typename T>
std::ostream& operator<<(std::ostream& os, const std::list<T>& lv)
{
    os << "(";
    //copy (lv.begin(), lv.end(), ostream_iterator<T>(os, ", "));
    //os << ")\n";
    
    for ( auto it= std::begin(lv), e= std::end(lv); it != e; ++it) {
        if ( distance(it,e) > distance(prev(e,1), e) )
            os <<  *it << " , ";
        else
            os << *it << ")\n";      
    }
    
    return os;
}

template <>
std::ostream& operator<<(std::ostream& os, const std::list<std::vector<int> >& lv)
{
    for ( auto it= std::begin(lv), e= std::end(lv); it != e; ++it) {
        os << "{ [";
        for ( auto vit = std::begin(*it), ve= std::end(*it); vit != ve; ++vit ) {
            if ( vit < ve - 1)
                os <<  *vit << " , ";
            else
                os << *vit << " ] }\n";    
        }  
    }
    return os;
}




template<typename F, typename T>
class Condition
{
    public:
        Condition(F& f) : f(f) {}
        T operator() (const T& x)
        {
            return f(x);
        }

    private:
        F f;
};

template <typename T>
T line(T m, T b, T x)
{
    return m*x + b;
}

template <typename T>
T xx(T a, T b, T x)
{
    return a*x*x + b;
}

template <typename T>
class linear_f
{
    public:
        linear_f(T a) : a(a) {}

        const T operator() (const T& x)
        {
            return a*x;
        }
    private:
        T a;
};

template<typename T>
class sinax_f
{
    public:
        sinax_f(T a) : a(a) {}

        const T operator() (const T& x)
        {
            return sin(x);
        }
    private:
        T a;
};

template <typename T>
class constant_f
{
    public:
        constant_f(T a) : a(a) {}

        const T operator() (const T& x, const T& y)
        {
            return a;
        }
    private:
        T a;
};

// Numerical Analysis by Burden & Faires, Functions Chapter 12.1
class U0y
{
    public:
    const double operator() (const double& x, const double& y)
    {
        return 0;
    }
};

class Ux0
{
    public:  // No constructor, not needed anyway because not initializing any variables
    const double operator() (const double& x, const double& y)
    {
        return 0;
    }
};

class Uxb
{
    public:
    const double operator() (const double& x, const double& y)
    {
        return 200*x;
    }
};

class Udy
{
    public:
    const double operator() (const double& x, const double& y)
    {
        return 200*y;
    }
};


template <typename T>
class x_e_to_y
{
    public:
        x_e_to_y();  // default constructor

        const T operator() (const T& x, const T& y)
        {
            return x*std::exp(y);
        }
    //private:
    //    T a;
};

class poissonIdx
{
    using VI = std::vector<int>;
    using VD = std::vector<double>;
    using VVI = std::vector<std::vector<int> >;

    //friend class poisson_dif_eqn_f;

    public:
        poissonIdx(VI i, VI j, VD constants)
        {
            i2check = i;
            j2check = j;
            c2assign = constants;
        }

        VI index1(int i, int j)
        {
            VI values = {i, j};
            return values;
        }
        
        VI index2(int i, int j)
        {
            VI values = {i+1, j};
            return values;
        }

        VI index3(int i, int j)
        {
            VI values = {i-1, j};
            return values;
        }

        VI index4(int i, int j)
        {
            VI values = {i, j-1};
            return values;
        }

        VI index5(int i, int j)
        {
            VI values = {i, j+1};
            return values;
        }

        //template <typename T, typename rhs2D_f>
        //VVI indexMap(poisson_dif_eqn_f<T, rhs2D_f> pde_f, VVI idx2move2rhs) 
        VVI indexMap(int numCols, VVI idx2move2rhs)
        {
            VVI newIdx = idx2move2rhs;
            /*
            for ( int i= 0; i < idx2move2rhs.size(); i++ ) {
                auto found1 = find( (i2check).begin(), (i2check).end(), idx2move2rhs[i][0] );
                auto found2 = find( (j2check).begin(), (j2check).end(), idx2move2rhs[i][1] );

                if ( found1 != idx2move2rhs[i].end() ) {
                    //newIdx[i][0] = 
                }
            }
            */
           for (int row=0; row < newIdx.size(); row++ ) {
                newIdx[row][0] = numCols - idx2move2rhs[row][1];
                newIdx[row][1] = idx2move2rhs[row][0];
               
           }

            return newIdx;
        }
        
        template <typename T>
        std::vector<T> which_const(std::vector<T> constants, VI cidx)
        {
            std::vector<T> thisconst;
            if ( cidx.empty() ) {
                thisconst.push_back(1);
                return thisconst;
            }
            for ( int i= 0; i < cidx.size(); i++ ) {
                thisconst.push_back( -1 * constants[cidx[i]] );
            }
            return thisconst;
        }
        

        std::tuple<VI,VVI> which_idx_2move(int i, int j)
        // Future : Variadic function? Or is this okay? I think it's okay :) 
        {
            VI cidx;
            allIndices = { index1(i,j) , index2(i,j), index3(i,j), index4(i,j), index5(i,j)};
            auto resulti = std::find(begin(i2check), end(i2check), allIndices[0][0] );  // initialize here to put on stack before loop
            auto resultj = std::find(begin(j2check), end(j2check), allIndices[0][1] );
            VVI idx2move2rhs;
            //for (auto i_it= std::begin(allIndices), i_e= std::end(allIndices); i_it != i_e; ++i_it ) {
            for (int i_it= 0; i_it < allIndices.size(); i_it++ ) {
                //for (auto j_it= std::begin(*i_it), j_e= std::end(*i_it); j_it != j_e; advance(j_it,1) ) {
                resulti = std::find(std::begin(i2check), std::end(i2check), allIndices[i_it][0] );
                resultj = std::find(begin(j2check), end(j2check), allIndices[i_it][1] );

                if ( (resulti != end(i2check) ) || (resultj != end(j2check)) ) {
                    idx2move2rhs.push_back( allIndices[i_it] );
                    cidx.push_back(i_it);
                }
                //}
            }
            return std::make_tuple(cidx, idx2move2rhs);
        }

        bool is_boundary_terms(int i, int j)
        {
            bool result;
            allIndices = { index1(i,j) , index2(i,j), index3(i,j), index4(i,j), index5(i,j)};
            //cout << "\nallIndices.size() = "<< allIndices.size() << ", allIndices: \n" << allIndices << endl;

            auto resulti = std::find(begin(i2check), end(i2check), allIndices[0][0] );  // initialize here to put on stack before loop
            auto resultj = std::find(begin(j2check), end(j2check), allIndices[0][1] );

            for (int i_it= 0; i_it < allIndices.size(); i_it++ ) {
                resulti = std::find(std::begin(i2check), std::end(i2check), allIndices[i_it][0] );
                resultj = std::find(begin(j2check), end(j2check), allIndices[i_it][1] );

                if ( (resulti != end(i2check) ) || (resultj != end(j2check)) ) {
                    result = true;
                    return true;
                }
                else { result = false; }
            }
            return result;
        }
        

    private:
        VI i2check;
        VI j2check;
        VD c2assign;
        VVI allIndices;
};

template <typename T,typename rhs2D_f>
class poisson_dif_eqn_f
{
    public:
        using VI = std::vector<int>;
        using VVI = std::vector<std::vector<int> >;
        //using predicate = std::function<rhs2D_f>;
        poisson_dif_eqn_f( const rhs2D_f f, std::vector<T>& x, std::vector<T>& y, int n=4, int m=4, std::vector<int>& bidx={0}, T h=0.125, T k=0.125) : _f(f), x(x), y(y), _n(n), _m(m), _h(h), _k(k), bidx(bidx)
        {
            assert( h == k ); // h must equal k for setting up matrix later with mtl::laplacian due to the 4.
            constants = {2*( pow(h/k,2) + 1) , -1, -1, -pow(h/k,2) , -pow(h/k,2) };
            
            auto left = U0y(); 
            auto bottom = Ux0();
            auto top =  Uxb();
            auto right = Udy();
            auto W = build_W<U0y, Ux0, Uxb, Udy>(left, bottom, top, right);

            poissonIdx pidx( {0,_n},{0,_m}, constants );

            int len_b = (m-1) * (n-1);
            //rhs_b( len_b );
            //rhs_b(len_b);
            //rhs_b = 0.0;

            rhs_b = build_rhs_b< poissonIdx >(x, y, n, m, W, bidx, f, pidx );
            cout << "constants: " << constants << endl;
            cout << "W :\n" << W;
            cout << "rhs_b : "  << rhs_b << endl;
        }

        template <typename B1, typename B2, typename B3, typename B4>
        //dense2D<T> build_W(B1 f1, B2 f2, B3 f3, B4 f4)
        std::vector<std::vector<T> > build_W(B1 left, B2 bottom, B3 top, B4 right)
        {
            // assign values to W, this is boundary conditions only for mesh in box arrangement
            //dense2D<T> W((_n+1), (_m+1));
            std::vector<std::vector<T> > W(_n+1, std::vector<T> (_m+1, 0.0) );  // automatically filled with vectors of zeros
            cout << "Initializing W:\n" << "W.size() = " << W.size() << " , W[0].size() = " << W[0].size() <<"\n"<< W << endl;
            //W = 0.0;
            for ( int i= 0; i < W.size(); i++ ) {
                for ( int j=0 ; j < W[0].size(); j++) {
                    if ( i==0 )
                        W[i][j] = top( x[j], y[ _m - i]  ); // y[4]
                    if ( j== 0 )
                        W[i][j] = left( x[j], y[_m - i] ); // x[0]
                    if ( j == ( _m ) )
                        W[i][j] = right( x[ _n ], y[_m - i] ); // x[4]
                    if ( i == ( _n ) )
                        W[i][j] = bottom( x[j], y[ _m - i] ); // y[0]
                }
            }
            return W;
        }

        template <typename BoundType>
        dense_vector<T> build_rhs_b(std::vector<T>& x , std::vector<T>& y, int n, int m, std::vector<std::vector<T> >& W, std::vector<int>& bidx, rhs2D_f f, BoundType pidx)
        {
            // BT = poissonIdx
            //BT pidx({0,_n},{0,_m}, constants);
            std::vector<T> thisconst;
            VI cidx;
            VVI which_idx;
            VVI mapped_idx;
            // assign values to rhsb
            int l=0;
            int len_b = (m-1) * (n-1);
            //dense_vector<T> rhs_b( len_b , 0.0 );
            dense_vector<T, mtl::vec::parameters<tag::col_major>> temp_rhs( len_b, 0.0 );

            for ( int j= m-1; j >= 1; j-- ){
                for ( int i= 1; i < n; i++ ) {
                    l= i + (m - 1 - j)*(n-1) - 1;

                    auto result1 = find(begin(bidx), end(bidx), i);
                    auto result2 = find(begin(bidx), end(bidx), j);
                    
                    //cout << "pidx.is_boundary_terms(" << i << ", "<< j<< ") = " << pidx.is_boundary_terms(i,j) << " , l = "<< l << endl;
                    //cout << " x : \n" << x;
                    //cout << " y : \n" << y;

                    if ( pidx.is_boundary_terms(i,j) || result1 != end(bidx) || result2 != end(bidx) ) {                          
                        auto [ cidx, which_idx ] = pidx.which_idx_2move(i,j);  // tuple assignment
                        thisconst = pidx.which_const(constants, cidx);
                        cout << "\nwhich_idx = \n" << which_idx << "\n";

                        mapped_idx = pidx.indexMap(_m, which_idx);
                        cout << "\nmapped_idx = \n" << mapped_idx << "\n";                        
                        
                        for ( int k= 0; k < mapped_idx.size(); k++ ) {
                            temp_rhs[l] +=  thisconst[k] * W[ mapped_idx[k][0] ][ mapped_idx[k][1] ] + (-pow(_h,2)) * f( x[i], y[j] ) ;
                            cout << "W[" << mapped_idx[k][0]<< "][" << mapped_idx[k][1] << "] =" << W[mapped_idx[k][0]][mapped_idx[k][1]] << " , f(" << x[i] << ", " << y[j] << ") = " <<  f(x[i],y[j]) << endl;
                            cout << "l = " << l << " ,   i = " << i << ", j = "<< j <<  endl;
                        }
                    }
                    else {
                        temp_rhs[l] += -pow(_h,2) *f( x[i], y[j] );
                    }
                }   
            }
            cout << "num_rows(temp_rhs) = " << num_rows(temp_rhs) << endl;
            cout << "num_cols(temp_rhs) = " << num_cols(temp_rhs) << endl;
            cout << "size(temp_rhs) = " << size(temp_rhs) << endl;

            return temp_rhs;
        }

        int getCols() {
            return(_m);
        }

        dense_vector<T, mtl::vec::parameters<tag::col_major>> get_rhsb() {
            return(rhs_b);
        }

        std::vector<T> get_const()
        {
            return constants;
        }
        

    private:
        T _h, _k;
        std::vector<T> x, y;
        dense_vector<T, mtl::vec::parameters<tag::col_major>> rhs_b;
        VI bidx;
        rhs2D_f _f;
        int _n, _m;
        std::vector<T> constants;
        //dense2D<T> W;
        std::vector<std::vector<T> > W;
};


template <typename Points, typename T>
class Boundary
{
    public:
    int size;
    mtl::dense_vector<Points> points;
    std::vector<std::function< T(T)> > functors;

    //Boundary(std::initializer_list<Points> plist, std::initializer_list<std::function<T(T)> > flist) : size(plist.size())
    Boundary(Points plist, std::function<T(T)> flist) : size(plist.size())
    {
        assert( size(plist) == size(flist) );
        points = plist;
        functors = flist;
    }

};

const double pi= 3.14159265;

int main(int argc, char* argv[])
{
    {    
    //double A0[][2]= {{3, 4}, {5, 6}};
    //const mtl::dense2D<double> A(A0);
    const mtl::dense2D<double> A = {{3, 4}, {5, 6}};
    cout << "A = " << A << endl;
    }

    cout << "\n";
    calc_l(4,4);
    cout << "\n";
    calc_lr(4,4);

    Point2D<double> pt(0.0, 2.0);
    cout << pt;

    PointND<double> pt2({2.7,100.5});
    cout << pt2;

    cout << PointND<double> {0.9999875};

    cout << PointND<double> {3.2222, 6.5};

    cout << PointND<double> {3.2222, 6.5, 55.555};

    cout << PointND<double> {100.1, 200.2, 300.3, 444.4};

    //mtl::dense_vector<PointND<double>> myPoints = {{0.0,0.5} ,{-3.14, 3.14} };
    //std::vector<PointsND<double>> myPoints = { {0.0,0.5} ,{-3.14, 3.14} };


    //Boundary<Condition<linear_f<double>,double>, PointND<double>, double> left({PointND<double>{0.0,0.5}} , {Condition<linear_f<double>,double>(linear_f<double>(1.0,2.2)) } );

    cout << "\n linear_f functor test\n";
    linear_f y(1.0);
    cout << y(0) << ", " << y(1.0) << ", " << y(333.567) << "\n\n";

    cout << "\n sinax functor test\n";
    sinax_f myfunc(1.0);
    cout << myfunc(pi/180*30) << ", " << myfunc(1.0) << ", " << myfunc(333.567) << "\n\n";


    cout << "\n\n Vector Functor test:\n---------------\n";

    std::vector<std::function<double()>> functors;
    //functors.push_back([&] { return 100.0; });   // works
    //functors.push_back([&] { return  20.0; });   // works 
    // functors = { [&] { return 100; }, [&] { return 20; }  }; // works

    using std::bind;
    functors.push_back( bind(line<double>,1.0,0.0,3.14) );   // works
    functors.push_back( bind(  xx<double>,1.0,0.0,2.0 ) );
    functors.push_back( bind( y, 3.14 ) );
    functors.push_back( bind( myfunc, 3.14/2 ) );
    functors.push_back( bind( linear_f(1.0), 3.14 ) );
    functors.push_back( bind( sinax_f(1.0), 3.14/2 ) );

    for ( auto f : functors )
    {
        cout << f() << endl;
    }

    cout << "\n\ninit_functor test:\n";
    std::vector<std::function<double()>> init_functors;
    init_functors = { bind(sinax_f(1.0), 3.14/2), bind(linear_f(1.0), 2.0 ) };

    for ( auto f : init_functors )
    {
        cout << f() << endl;
    }

    using namespace std::placeholders;
    cout << "\n\nplace_holder functor test:\n-------------\n";

    auto pfun = bind(sinax_f<double>(1.0), _1);
    double x0 = 0.0;
    double h  = 0.1;
    for (int i=0; i < 40; i++)
    {
        cout << "pfun(" << x0+i*h << ") = " << pfun(x0+i*h) << endl;
    }

    
    cout << "\n\nVector place_holder test:\n-------------\n";
    std::vector<std::function< double(double) >> ph_functors;
    //ph_functors = { bind(sinax_f<double>(1.0), _1), bind(linear_f<double>(1.0), _1 ) };  // works
    //ph_functors = { bind(y, _1), bind( myfunc, _1 ) };   // also works
    ph_functors = { y , myfunc };  // works again!

    for ( auto phf : ph_functors )
    {
        cout << phf(3.14/2) << endl;
    }


    
    cout << "\n\nmtl_functor test 1:\n---------------\n";
    dense_vector<double> xx;
    xx = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0};

    dense_vector<std::function<double(double)>> mtl_functor(2);
    mtl_functor[0] = y; //{ y, myfunc };
    mtl_functor[1] = myfunc;

    for ( auto f : mtl_functor )
    {
        cout << f( pi/2 ) << endl;
    }

    cout << "\n\nmtl_functor test 2:\n---------------\n";
    dense_vector<std::function<double(double)>> mtl_init_f({ y , myfunc });
    /* 
    mtl_init_f = { y , myfunc };  
    */

    for ( auto f : mtl_init_f )
    {
        for (int i=0; i < 40; i++) {
            cout << f( x0+i*h ) << endl;
        }
        
    }
    
    



    cout << "\n\n---------------------\nBoundary Test:\n--------------------\n";
    
    mtl::dense_vector<double , mtl::vec::parameters<mtl::tag::col_major> > x(5);
    x = 5.0;

    
    //cout << "size of mtl vector is: " << mtl::traits::num_cols<mtl::dense_vector<double>>(x) << endl;
    //cout << "size of mtl vector is: " << mtl::traits::num_rows(x) << endl;
    cout << "size of mtl vector is: " << mtl::size(x) << endl;   // calls mtl::traits::size() constructors
    
    //Boundary left( PointND<double>{0.0,0.5}, y );
    using namespace mtl;
    
    cout << "size of mtl vector is: " << num_cols(x) << endl;
    cout << "size of mtl vector is: " << num_rows(x) << endl;
    cout << "size of mtl vector is: " << size(x) << endl;


    cout << "\n\nPoisson Solver:\n----------------\n";
    {
        const int n= 4, m= 4;
        double xa= 0, xb= 0.5, yc= 0, yd= 0.5;
        double h = (xb-xa)/n, k= (yd-yc)/m;
        //dense_vector<double> x(n+1, 0.0), y(m+1, 0.0);
        std::vector<double> x(n+1), y(m+1);
        for ( int i=0; i <= n; i++) {
            x[i] = xa + i*h;
        }
        cout << "x is :" << x << endl;
        
        for ( int j=0; j <= m; j++) {
            y[j] = yc + j*k;
        }
        cout << "y is :" << y << endl;

        compressed2D<double>   A((n-1)*(m-1), (n-1)*(m-1));
        // Laplace operator discretized on a 3x3 grid
        mat::laplacian_setup(A, n-1, m-1);
        cout << "A is \n" << A;
        cout << "size(A) is : " << mtl::mat::size(A) << endl;  // or just size(A)
        cout << "num_cols(A) is : " << mtl::mat::num_cols(A) << endl; 
        cout << "num_rows(A) is : " << num_rows(A) << endl;

        dense_vector<double> b((n-1)*(m-1) , 0.0);
        int l= 0;
        for ( int j= m-1; j >= 1; j-- ){
            for ( int i= 1; i < n; i++ ) {
                l= i + (m - 1 - j)*(n-1) - 1;
                b[l] = constant_f( 0.0 )(x[i] ,y[j]);  // works
                // b[l] = [alpha=2.2222](double x, double y){ return alpha; }(x[i] ,y[j]); // works
            }
        }
        cout << "b is :" << b << "\n\n";

        /*  Store indices that are part of boundary, or are known, instead of Boundary class */

        // single list
        cout << "single list:\n";
        std::list<int> list1 = {0, 1, 2, 3, 4, 5};
        for (std::list<int>::iterator it= begin(list1), e= end(list1); it != e; ++it) {
            cout << *it << endl;
        }

        cout << "single list 2:\n";
        cout << list1;


        // list of lists
        cout << "\nlist of lists:\n";
        std::list<std::list<int> > list2 = {{0,1}, {0,2}, {0,3}, {1,0} };
        for (std::list<std::list<int> >::iterator it1= begin(list2), e1= end(list2); it1 != e1; ++it1) {
            for(std::list<int>::iterator it2= begin(*it1), e2= end(*it1); it2 != e2; std::advance(it2, 2)  ) {
                cout << "{ "<< *it2 << "," <<  *std::next(it2,1) << " }"  << endl;
            }
        }

        cout << "\nlist of vectors of known size:\n";
        using std::list;
        using std::vector;
        using std::advance;
        list<vector<int> > lv;
        lv.push_back(vector<int>{1,2});
        lv.push_back(vector<int>{3,3});
        lv.push_back(vector<int>{9,8});

        for (list<vector<int> >::iterator it= begin(lv), e= end(lv); it != e; ++it) {
            
            for ( auto vit = begin(*it), ve= end(*it); vit != ve; advance(vit,2) ) {
                cout << "( " << *vit << " , " << *next(vit,1) << " )\n";
            }
            
        }

        cout << "\n\nPrint a List of vectors:---------------------\n";
        cout << lv;

        cout << "\n\nvectors of vectors of known size\n";
        vector< vector<int> > vv;
        vv = {{0,0},{0,1},{0,2}};
        vv.push_back({1,0});
        
        for ( auto it= std::begin(vv), e= std::end(vv); it < e; ++it ) {
            for ( auto vi= std::begin(*it), ve= std::end(*it); vi < ve; ++vi ) {
                if ( vi < ve - 1)
                    cout << "[" << *vi << " , ";
                else
                    cout << *vi << "]\n";
            }
        }

        cout << "\n\nVectors of vectors 2:\n";
        for ( auto& it: vv) {
            cout << "[";
            for ( auto& vi: it) {
                cout << vi << " , ";
            }
            cout << "]" << endl;
        }

        cout << "\n\nVector of Vectors 3" << endl;
        cout << vv;

        cout << "\n\nVector of ints" << endl;
        std::vector<int> v1= {1,2,3,4,5,6,7,8,9};
        cout << v1;

        // Build vector b in Ax = b
        // Function that accepts 1) difference equation Functor 2) list of boundary point indices, 
        // and known values for inner boundaries. Returns vector b.
        
        dense2D<double> W((n+1),(m+1));
        W = 0.0;
        cout << "\n\n" << W;

        cout << "\n\nconstant_f :";
        constant_f<double> One(1.0);
        cout << One(0,0);

        cout << "\n\nPoisson Test:\n";
        constant_f zeroFunc(0.0);
        vector<int> bidx = {0};
        poisson_dif_eqn_f poissonSolver( zeroFunc, x, y, n, m, bidx, h, k );
        dense_vector<double, mtl::vec::parameters<tag::col_major>> rhs_b = poissonSolver.get_rhsb();
        cout << "poissonSolver.get_rhsb() = " << poissonSolver.get_rhsb() << endl;

        /* From SimuNova Tutorial for MTL4  */
        using namespace itl;
        typedef compressed2D<double>  matrix_type;
        // Create an ILU(0) preconditioner
        pc::ilu_0<matrix_type> P(A);

        dense_vector<double, mtl::vec::parameters<tag::col_major>> solution(num_rows(rhs_b) , 0.0);
        cout << "num_rows(rhs_b) = " << num_rows(rhs_b) << endl;
        cout << "size(A) = " << size(A) << endl;
        cout << "get_const() = " << poissonSolver.get_const();

        // Termination criterion: r < 1e-6 * b or N iterations
        noisy_iteration<double>       iter(b, 500, 1.e-6);

        // Solve Ax == b with left preconditioner P
        bicgstab(A, solution, rhs_b, P, iter);

        cout << solution << endl;




    }// end of Poisson Equation Solver
        



    return 0;
}
