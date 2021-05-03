#include <Rcpp.h>
using namespace Rcpp;

// Function from StackOverflow - 
// From SO https://stackoverflow.com/questions/59907035/memory-efficient-method-to-create-dist-object-from-distance-matrix
// Users @dww & @Alexis 

// [[Rcpp::export]]
NumericVector as_dist(const NumericMatrix& mat){
#include <cstddef> // size_t
    
    std::size_t nrow = mat.nrow();
    std::size_t ncol = mat.ncol();
    std::size_t size = nrow * (nrow - 1) / 2;
    NumericVector ans(size);
    
    if (nrow > 1) {
        std::size_t k = 0;
        for (std::size_t j = 0; j < ncol; j++) {
            for (std::size_t i = j + 1; i < nrow; i++) {
                ans[k++] = mat(i,j);
            }
        }
    }
    
    ans.attr("class") = "dist";
    ans.attr("Size") = nrow;
    ans.attr("Diag") = false;
    ans.attr("Upper") = false;
    return ans;
}