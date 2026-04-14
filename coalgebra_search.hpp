#pragma once

#include <optional>
#include <unordered_map>
#include <unordered_set>

#include "DegreeZeroGraphWedge.hpp"
#include "GC.hpp"

namespace coalgebra_search {

template <typename LeftBasis, typename RightBasis>
class DirectSumBasis {
	public:
		enum class Tag : unsigned char { Left, Right };

		DirectSumBasis() = default;

		static DirectSumBasis from_left(const LeftBasis& value) {
			DirectSumBasis basis;
			basis.tag_ = Tag::Left;
			basis.left_ = value;
			return basis;
		}

		static DirectSumBasis from_right(const RightBasis& value) {
			DirectSumBasis basis;
			basis.tag_ = Tag::Right;
			basis.right_ = value;
			return basis;
		}

		Tag tag() const noexcept { return tag_; }
		const LeftBasis& left() const noexcept { return left_; }
		const RightBasis& right() const noexcept { return right_; }

		signedInt compare(const DirectSumBasis& other) const noexcept {
			if (tag_ != other.tag_) {
				return static_cast<signedInt>(tag_) < static_cast<signedInt>(other.tag_) ? -1 : 1;
			}

			if (tag_ == Tag::Left) {
				return left_.compare(other.left_);
			}
			return right_.compare(other.right_);
		}

		bool operator<(const DirectSumBasis& other) const noexcept {
			return compare(other) < 0;
		}

		bool operator==(const DirectSumBasis& other) const noexcept {
			if (tag_ != other.tag_) {
				return false;
			}
			if (tag_ == Tag::Left) {
				return left_ == other.left_;
			}
			return right_ == other.right_;
		}

	private:
		Tag tag_ = Tag::Left;
		LeftBasis left_{};
		RightBasis right_{};
};

template <typename GraphGC, typename LeftGC, typename RightGC>
using TwoFactorWedge = DegreeZeroGraphWedge<
	static_cast<Int>(2 * GraphGC::GraphType::N_EDGES_),
	2
>;

template <typename GraphGC>
using OneFactorWedge = DegreeZeroGraphWedge<
	static_cast<Int>(2 * GraphGC::GraphType::N_EDGES_),
	1
>;

template <typename GraphGC, typename LeftGC, typename RightGC>
using JointConstraintBasis = DirectSumBasis<
	typename GraphGC::ContGraphType,
	TwoFactorWedge<GraphGC, LeftGC, RightGC>
>;

template <typename LeftGC, typename RightGC, typename WedgeType>
typename WedgeType::LinComb wedge_target(const LeftGC& left, const RightGC& right) {
	using WedgeLinComb = typename WedgeType::LinComb;

	WedgeLinComb target;
	for (const auto& left_be : left.data()) {
		for (const auto& right_be : right.data()) {
			auto wedge = WedgeType::from_graphs(left_be.getValue(), right_be.getValue());
			target.append_in_basis_order(
				wedge,
				left_be.getCoefficient() * right_be.getCoefficient()
			);
		}
	}

	target.standardize_and_sort();
	return target;
}

template <typename GraphGC, typename LeftGC, typename RightGC>
typename TwoFactorWedge<GraphGC, LeftGC, RightGC>::LinComb cobracket_of_gc(const GraphGC& gc) {
	using OneWedgeType = OneFactorWedge<GraphGC>;
	using WedgeType = TwoFactorWedge<GraphGC, LeftGC, RightGC>;
	using WedgeLinComb = typename WedgeType::LinComb;

	WedgeLinComb result;
	for (const auto& be : gc.data()) {
		auto graph_wedge = OneWedgeType::from_graph(be.getValue()).cobracket();
		graph_wedge.scalar_multiply(be.getCoefficient());
		result += graph_wedge;
	}

	result.standardize_and_sort();
	return result;
}

template <typename GraphGC, typename LeftGC, typename RightGC>
std::unordered_map<
	typename GraphGC::GraphType,
	VectorSpace::LinComb<JointConstraintBasis<GraphGC, LeftGC, RightGC>, fieldType>
> build_joint_constraint_map(
	const std::unordered_set<typename GraphGC::GraphType>& support,
	bool verbose = false
) {
	using GraphType = typename GraphGC::GraphType;
	using ContGraphType = typename GraphGC::ContGraphType;
	using OneWedgeType = OneFactorWedge<GraphGC>;
	using WedgeType = TwoFactorWedge<GraphGC, LeftGC, RightGC>;
	using JointBasis = JointConstraintBasis<GraphGC, LeftGC, RightGC>;
	using JointLinComb = VectorSpace::LinComb<JointBasis, fieldType>;

	std::unordered_map<GraphType, JointLinComb> map;
	map.reserve(support.size());

	std::size_t processed = 0;
	for (const auto& graph : support) {
		JointLinComb image;

		auto differential = GraphGC(graph, AssumeBasisOrderTag{}).d_contraction();
		differential.standardize_all();
		differential.sort_elements();
		for (const auto& be : differential.data()) {
			image.append_in_basis_order(
				JointBasis::from_left(be.getValue()),
				be.getCoefficient()
			);
		}

		auto cobracket = OneWedgeType::from_graph(graph).cobracket();
		for (const auto& be : cobracket) {
			image.append_in_basis_order(
				JointBasis::from_right(be.getValue()),
				be.getCoefficient()
			);
		}

		image.sort_elements();
		map.emplace(graph, std::move(image));

		++processed;
		if (verbose && processed % 50 == 0) {
			std::cout << "built joint rows: " << processed
				  << " / " << support.size() << std::endl;
		}
	}

	return map;
}

template <typename GraphGC, typename LeftGC, typename RightGC>
std::optional<GraphGC> find_cobracket_lift_in_support(
	const std::unordered_set<typename GraphGC::GraphType>& support,
	const LeftGC& left,
	const RightGC& right,
	bool verbose = false
) {
	using GraphType = typename GraphGC::GraphType;
	using ContGraphType = typename GraphGC::ContGraphType;
	using WedgeType = TwoFactorWedge<GraphGC, LeftGC, RightGC>;
	using JointBasis = JointConstraintBasis<GraphGC, LeftGC, RightGC>;
	using JointLinComb = VectorSpace::LinComb<JointBasis, fieldType>;

	static_assert(
		static_cast<Int>(2 * LeftGC::GraphType::N_EDGES_ + 2 * RightGC::GraphType::N_EDGES_)
			== static_cast<Int>(2 * GraphGC::GraphType::N_EDGES_),
		"wedge target must match the total half-edge count of the search space"
	);

	const auto target_wedge = wedge_target<LeftGC, RightGC, WedgeType>(left, right);
	JointLinComb target;
	for (const auto& be : target_wedge) {
		target.append_in_basis_order(
			JointBasis::from_right(be.getValue()),
			be.getCoefficient()
		);
	}
	target.sort_elements();

	if (verbose) {
		std::cout << "target wedge terms: " << target_wedge.size() << std::endl;
	}

	const auto joint_map = build_joint_constraint_map<GraphGC, LeftGC, RightGC>(support, verbose);
	if (verbose) {
		std::cout << "built joint constraint map with " << joint_map.size() << " domain graphs" << std::endl;
	}

	VectorSpace::wiedemann_primitive_finder<JointBasis, GraphType, fieldType> solver(joint_map);
	auto solution_opt = solver.find_primitive_or_empty(target);
	if (!solution_opt.has_value()) {
		return std::nullopt;
	}

	auto solution = *solution_opt;
	solution.standardize_and_sort();
	GraphGC candidate(std::move(solution.raw_elements_nonconst()));

	auto d_candidate = candidate.d_contraction();
	d_candidate.standardize_all();
	d_candidate.sort_elements();
	if (d_candidate.data().size() != 0) {
		return std::nullopt;
	}

	const auto candidate_wedge = cobracket_of_gc<GraphGC, LeftGC, RightGC>(candidate);
	if (!(candidate_wedge == target_wedge)) {
		return std::nullopt;
	}

	return candidate;
}

} // namespace coalgebra_search

namespace std {

template <typename LeftBasis, typename RightBasis>
struct hash<coalgebra_search::DirectSumBasis<LeftBasis, RightBasis>> {
	std::size_t operator()(
		const coalgebra_search::DirectSumBasis<LeftBasis, RightBasis>& basis
	) const noexcept {
		std::size_t seed = std::hash<unsigned char>{}(static_cast<unsigned char>(basis.tag()));
		const auto mix = [&seed](std::size_t value) {
			seed ^= value + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2);
		};

		if (basis.tag() == coalgebra_search::DirectSumBasis<LeftBasis, RightBasis>::Tag::Left) {
			mix(std::hash<LeftBasis>{}(basis.left()));
		} else {
			mix(std::hash<RightBasis>{}(basis.right()));
		}

		return seed;
	}
};

} // namespace std
