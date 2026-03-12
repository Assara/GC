#pragma once

#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "VectorSpace/LinComb.hpp"
#include "graph.hpp"

class GraphVisualizer {
	public:
		struct Options {
			std::string graph_name = "G";
			bool show_edge_labels = true;
			bool show_valence_labels = true;
			bool color_by_valence = true;
		};

		template <
			Int N_VERTICES,
			Int N_EDGES,
			Int N_OUT_HAIR,
			Int N_IN_HAIR,
			signedInt c,
			signedInt d,
			typename fieldType
		>
		static std::string to_dot(
			const Graph<N_VERTICES, N_EDGES, N_OUT_HAIR, N_IN_HAIR, c, d, fieldType>& graph,
			const Options& options
		) {
			std::ostringstream out;
			const auto valences = graph.valence_array();

			out << "graph " << options.graph_name << " {\n";
			out << "  layout=neato;\n";
			out << "  overlap=false;\n";
			out << "  splines=true;\n";
			out << "  node [shape=circle, style=filled, fillcolor=white, fontname=\"Helvetica\"];\n";
			out << "  edge [fontname=\"Helvetica\"];\n";

			for (Int v = 0; v < N_VERTICES; ++v) {
				out << "  v" << +v << " [label=\"";
				if (options.show_valence_labels) {
					out << +v << " (" << +valences[v] << ")";
				} else {
					out << +v;
				}
				out << "\"";

				if (options.color_by_valence) {
					out << ", fillcolor=\"" << color_for_valence(valences[v]) << "\"";
				}

				out << "];\n";
			}

			for (Int e = 0; e < N_EDGES; ++e) {
				auto [u, v] = graph.getEdge(e);
				out << "  v" << +u << " -- v" << +v;
				if (options.show_edge_labels) {
					out << " [label=\"" << +e << "\"]";
				}
				out << ";\n";
			}

			out << "}\n";
			return out.str();
		}

		template <
			Int N_VERTICES,
			Int N_EDGES,
			Int N_OUT_HAIR,
			Int N_IN_HAIR,
			signedInt c,
			signedInt d,
			typename fieldType
		>
		static std::string to_dot(
			const Graph<N_VERTICES, N_EDGES, N_OUT_HAIR, N_IN_HAIR, c, d, fieldType>& graph
		) {
			return to_dot(graph, Options{});
		}

		template <
			Int N_VERTICES,
			Int N_EDGES,
			Int N_OUT_HAIR,
			Int N_IN_HAIR,
			signedInt c,
			signedInt d,
			typename fieldType
		>
		static void write_dot(
			const std::filesystem::path& output_path,
			const Graph<N_VERTICES, N_EDGES, N_OUT_HAIR, N_IN_HAIR, c, d, fieldType>& graph,
			const Options& options
		) {
			std::ofstream out(output_path);
			out << to_dot(graph, options);
		}

		template <
			Int N_VERTICES,
			Int N_EDGES,
			Int N_OUT_HAIR,
			Int N_IN_HAIR,
			signedInt c,
			signedInt d,
			typename fieldType
		>
		static void write_dot(
			const std::filesystem::path& output_path,
			const Graph<N_VERTICES, N_EDGES, N_OUT_HAIR, N_IN_HAIR, c, d, fieldType>& graph
		) {
			write_dot(output_path, graph, Options{});
		}

		template <typename GraphType, typename fieldType>
		static void write_lincomb_dot_files(
			const std::filesystem::path& output_dir,
			const VectorSpace::LinComb<GraphType, fieldType>& lincomb,
			const std::string& prefix,
			const Options& options
		) {
			std::filesystem::create_directories(output_dir);

			std::size_t index = 0;
			for (const auto& be : lincomb) {
				Options current = options;
				current.graph_name = prefix + "_" + std::to_string(index);

				const auto output_path = output_dir / (current.graph_name + ".dot");
				std::ofstream out(output_path);
				out << "// coefficient: " << be.getCoefficient() << "\n";
				out << to_dot(be.getValue(), current);
				++index;
			}
		}

		template <typename GraphType, typename fieldType>
		static void write_lincomb_dot_files(
			const std::filesystem::path& output_dir,
			const VectorSpace::LinComb<GraphType, fieldType>& lincomb,
			const std::string& prefix = "graph"
		) {
			write_lincomb_dot_files(output_dir, lincomb, prefix, Options{});
		}

	private:
		static const char* color_for_valence(Int valence) {
			switch (valence) {
				case 0: return "#f5f5f5";
				case 1: return "#d9edf7";
				case 2: return "#dff0d8";
				case 3: return "#fcf8e3";
				case 4: return "#f2dede";
				default: return "#e8d5ff";
			}
		}
};
