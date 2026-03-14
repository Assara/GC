
#pragma once

#include "LinComb.hpp"
#include "wiedemann_helper.hpp"



#include <vector>
#include <unordered_map>
#include <optional>
#include <memory>
#include <algorithm>
#include <iostream>

namespace VectorSpace {

	template<typename A, typename B, typename k>
		class wiedemann_primitive_finder {
			private:
				std::vector<B> domain_space_enumeration;               // row index -> B basis element  
				lil_matrix<k> map_representative;
				// map B -> A
				std::unordered_map<A, std::size_t> image_space_enumeration; // A basis -> column index


				using DenseDomainVec = std::unique_ptr<k[]>; //Assumed to have the size of the map_representative matrix
				using DenseImageVec = std::unique_ptr<k[]>; 


				std::size_t image_space_to_number(const A& a) {
					auto [it, inserted] =
						image_space_enumeration.try_emplace(a, image_space_enumeration.size());

					return it->second;
				}

				LinComb<size_t, k> map_to_enumeration_basis(const LinComb<A,k>& input) {
					LinComb<size_t, k> result;

					for (const auto& be : input ) {
						result.append_in_basis_order(BasisElement(image_space_to_number(be.getValue()), be.getCoefficient()));
					}


					result.sort_without_deduplicate();
					return result;
				}

				template<typename Pred>
				LinComb<size_t, k> map_to_enumeration_basis_filtered(const LinComb<A,k>& input, const Pred& predicate) {
					LinComb<size_t, k> result;

					for (const auto& be : input ) {
						if (!predicate(be.getValue())) {
							continue;
						}
						result.append_in_basis_order(BasisElement(image_space_to_number(be.getValue()), be.getCoefficient()));
					}

					result.sort_without_deduplicate();
					return result;
				}


				DenseImageVec map_to_enumeration_basis_dense(const LinComb<A,k>& input) {
					DenseImageVec result = map_representative.make_dense_image_vec_zero();

					for (const auto& be : input ) {
						auto it = image_space_enumeration.find(be.getValue());
						if (it == image_space_enumeration.end()) {
							// Project away terms outside the currently enumerated image basis.
							continue;
						}
						result[it->second] += be.getCoefficient();	

					}
					return result;
				}


				LinComb<B, k> domain_enumeration_inverse_from_dense(const DenseDomainVec& input) {
					LinComb<B, k> result;


					for (size_t i = 0; i< map_representative.domain_dim(); i++) {
						if (input[i] == k{0}) continue;
						result.append_in_basis_order(BasisElement(domain_space_enumeration[i], input[i]));
					}

					result.sort_without_deduplicate(); //this does not have to be sorted for current applications
					return result;
				}

				template<typename DeltaRange, typename ColumnBuilder>
				void build_from_delta_range(const DeltaRange& deltas, ColumnBuilder&& build_column) {
					std::size_t total_domain_size = 0;
					for (const auto& delta : deltas) {
						total_domain_size += delta.size();
					}
					domain_space_enumeration.reserve(total_domain_size);
					size_t n_matrix_entries = 0;

					for (const auto& delta : deltas) {
						for (const auto& entry : delta) {
							auto col = build_column(entry.second);
							domain_space_enumeration.push_back(entry.first);
							n_matrix_entries += col.size();
							map_representative.add_col(std::move(col));
						}
					}

					cout << "created sparse solver: domain_dim = " << domain_space_enumeration.size() << endl
						<< "image_space_dim = " <<  image_space_enumeration.size() << endl
						<< "number of matrix entries =" <<  n_matrix_entries << endl;
				}

				template<typename DeltaRange, typename ColumnBuilder>
				void build_from_delta_range_filtered(
					const DeltaRange& deltas,
					ColumnBuilder&& build_column,
					std::optional<std::size_t> max_non_zero_items = std::nullopt
				) {
					std::size_t total_domain_size = 0;
					for (const auto& delta : deltas) {
						total_domain_size += delta.size();
					}
					domain_space_enumeration.reserve(total_domain_size);
					size_t n_matrix_entries = 0;
					size_t skipped_dense_columns = 0;

					for (const auto& delta : deltas) {
						for (const auto& entry : delta) {
							auto col = build_column(entry.second);
							if (max_non_zero_items.has_value() && col.size() > *max_non_zero_items) {
								++skipped_dense_columns;
								continue;
							}

							domain_space_enumeration.push_back(entry.first);
							n_matrix_entries += col.size();
							map_representative.add_col(std::move(col));
						}
					}

					cout << "created filtered sparse solver: domain_dim = " << domain_space_enumeration.size() << endl
						<< "image_space_dim = " <<  image_space_enumeration.size() << endl
						<< "number of matrix entries =" <<  n_matrix_entries << endl;
					if (max_non_zero_items.has_value()) {
						cout << "skipped dense columns = " << skipped_dense_columns << endl;
					}
				}


			public:

				// -----------------------------------------------------
				// Constructor: build sparse matrix for δ : B → LinComb<A,k>
				// -----------------------------------------------------
				wiedemann_primitive_finder(const std::unordered_map<B, LinComb<A,k>>& delta) {
					domain_space_enumeration.reserve(delta.size());
					size_t n_matrix_entries = 0;

					for (const auto& entry : delta) {
						domain_space_enumeration.push_back(entry.first);

						auto col = map_to_enumeration_basis(entry.second);

						n_matrix_entries += col.size();
						map_representative.add_col(std::move(col));
					}

					cout << "created sparse solver: domain_dim = " << domain_space_enumeration.size() << endl
						<< "image_space_dim = " <<  image_space_enumeration.size() << endl
						<< "number of matrix entries =" <<  n_matrix_entries << endl;
				}

				wiedemann_primitive_finder(const std::vector<std::unordered_map<B, LinComb<A,k>>>& deltas) {
					build_from_delta_range(
						deltas,
						[this](const LinComb<A,k>& input) { return map_to_enumeration_basis(input); }
					);
				}

				template<typename Pred>
				static wiedemann_primitive_finder create_filtered(const std::unordered_map<B, LinComb<A,k>>& delta, const Pred& predicate) {
					return wiedemann_primitive_finder(delta, predicate);
				}

				template<typename Pred>
				static wiedemann_primitive_finder create_filtered(
					const std::vector<std::unordered_map<B, LinComb<A,k>>>& deltas,
					const Pred& predicate
				) {
					return wiedemann_primitive_finder(deltas, predicate);
				}

				template<typename Pred>
				static wiedemann_primitive_finder create_filtered(
					const std::unordered_map<B, LinComb<A,k>>& delta,
					const Pred& predicate,
					std::size_t max_non_zero_items
				) {
					return wiedemann_primitive_finder(delta, predicate, max_non_zero_items);
				}

				template<typename Pred>
				static wiedemann_primitive_finder create_filtered(
					const std::vector<std::unordered_map<B, LinComb<A,k>>>& deltas,
					const Pred& predicate,
					std::size_t max_non_zero_items
				) {
					return wiedemann_primitive_finder(deltas, predicate, max_non_zero_items);
				}

				template<typename Pred>
				wiedemann_primitive_finder(const std::unordered_map<B, LinComb<A,k>>& delta, const Pred& predicate) {
					domain_space_enumeration.reserve(delta.size());
					size_t n_matrix_entries = 0;

					for (const auto& entry : delta) {
						domain_space_enumeration.push_back(entry.first);

						auto col = map_to_enumeration_basis_filtered(entry.second, predicate);
						n_matrix_entries += col.size();
						map_representative.add_col(std::move(col));
					}

					cout << "created filtered sparse solver: domain_dim = " << domain_space_enumeration.size() << endl
						<< "image_space_dim = " <<  image_space_enumeration.size() << endl
						<< "number of matrix entries =" <<  n_matrix_entries << endl;
				}

				template<typename Pred>
				wiedemann_primitive_finder(
					const std::vector<std::unordered_map<B, LinComb<A,k>>>& deltas,
					const Pred& predicate
				) {
					build_from_delta_range_filtered(
						deltas,
						[this, &predicate](const LinComb<A,k>& input) {
							return map_to_enumeration_basis_filtered(input, predicate);
						}
					);
				}

				template<typename Pred>
				wiedemann_primitive_finder(
					const std::unordered_map<B, LinComb<A,k>>& delta,
					const Pred& predicate,
					std::size_t max_non_zero_items
				) {
					domain_space_enumeration.reserve(delta.size());
					size_t n_matrix_entries = 0;
					size_t skipped_dense_columns = 0;

					for (const auto& entry : delta) {
						auto col = map_to_enumeration_basis_filtered(entry.second, predicate);
						if (col.size() > max_non_zero_items) {
							++skipped_dense_columns;
							continue;
						}

						domain_space_enumeration.push_back(entry.first);
						n_matrix_entries += col.size();
						map_representative.add_col(std::move(col));
					}

					cout << "created filtered sparse solver: domain_dim = " << domain_space_enumeration.size() << endl
						<< "image_space_dim = " <<  image_space_enumeration.size() << endl
						<< "number of matrix entries =" <<  n_matrix_entries << endl
						<< "skipped dense columns = " << skipped_dense_columns << endl;
				}

				template<typename Pred>
				wiedemann_primitive_finder(
					const std::vector<std::unordered_map<B, LinComb<A,k>>>& deltas,
					const Pred& predicate,
					std::size_t max_non_zero_items
				) {
					build_from_delta_range_filtered(
						deltas,
						[this, &predicate](const LinComb<A,k>& input) {
							return map_to_enumeration_basis_filtered(input, predicate);
						},
						max_non_zero_items
					);
				}


				//delta_0 and delta_1 must have non intersected domains.
				wiedemann_primitive_finder(const std::unordered_map<B, LinComb<A,k>>& delta_0, const std::unordered_map<B, LinComb<A,k>>& delta_1) {
					domain_space_enumeration.reserve(delta_0.size() + delta_1.size());
					size_t n_matrix_entries = 0;

					for (const auto& entry : delta_0) {

						domain_space_enumeration.push_back(entry.first);   // each B becomes a row

						auto col = map_to_enumeration_basis(entry.second);

						n_matrix_entries += col.size();
						map_representative.add_col(std::move(col));	

					}

					for (const auto& entry : delta_1) {

						domain_space_enumeration.push_back(entry.first);   // each B becomes a row

						auto col = map_to_enumeration_basis(entry.second);

						n_matrix_entries += col.size();
						map_representative.add_col(std::move(col));	

					}

					cout << "created sparse solver: domain_dim = " << domain_space_enumeration.size() << endl
						<< "image_space_dim = " <<  image_space_enumeration.size() << endl
						<< "number of matrix entries =" <<  n_matrix_entries << endl;

				}


				std::optional<LinComb<B,k>> find_primitive_or_empty(LinComb<A,k> y) {
					cout << "using wierdemann_primitive_finder for LinComb:" << endl;

					auto y_enumerated = map_to_enumeration_basis_dense(y);

					cout << "mapped to enumeration!!!" << endl;

					//we do not need this anymore, but we do need heap space.
					image_space_enumeration.clear();
					wiedemann_solver<k> solver(map_representative);

					cout << "created solver!" <<  std::endl;

					std::optional<DenseImageVec> X_enumerated = solver.solve_MX_equals_y(y_enumerated);

					cout << "solved!!!" << endl;

					if (!X_enumerated.has_value()) {
						cout << "no solution!" << endl;

						return std::nullopt;
					}

					LinComb<B,k> X = domain_enumeration_inverse_from_dense(*X_enumerated);

					return std::optional(std::move(X));
				}

		};

} // namespace VectorSpace
