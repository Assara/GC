#include "GraphHomology/ContractionMatrices.hpp"
#include <iostream>
int main(){
 using namespace GraphHomology;
 const char* path="output/triangle_comparison_L10_cut_splits_run1/splits_L6";
 auto source=EnumeratedBasis<6,10,Parity::even>::load(path);
 auto target=EnumeratedBasis<6,9,Parity::even>::load(path);
 auto matrix=contraction_matrix(source,target);
 std::vector<int> sums(target.size()),parents(target.size());
 for(std::size_t c=0;c<matrix.columns();++c)for(auto i=matrix.offsets[c];i<matrix.offsets[c+1];++i){
  auto r=matrix.row_indices[i];sums[r]+=std::abs(int(matrix.coefficients[i]));++parents[r];
  if(std::abs(int(matrix.coefficients[i]))>=15){
   std::cout<<"source="<<c<<" target="<<r<<" coefficient="<<int(matrix.coefficients[i])<<" source_edges=";
   for(int e=0;e<15;++e){auto[a,b]=source.graphs[c].getEdge(e);std::cout<<int(a)<<','<<int(b)<<' ';}std::cout<<'\n';
  }
 }
 for(std::size_t r=0;r<target.size();++r){
  int splits=0;for(auto d:target.graphs[r].valence_array())splits+=(1<<(d-1))-d-1;
  if(sums[r]>splits){std::cout<<"row="<<r<<" splits="<<splits<<" parents="<<parents[r]<<" absolute_row_sum="<<sums[r]<<" degrees=";
   for(auto d:target.graphs[r].valence_array())std::cout<<int(d)<<',';std::cout<<'\n';}
 }
}
