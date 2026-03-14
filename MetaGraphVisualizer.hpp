#pragma once

#include <filesystem>
#include <fstream>
#include <functional>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>

#include "MetaGraph.hpp"

class MetaGraphVisualizer {
	public:
		struct Options {
			std::string graph_name = "MetaGraph";
			bool show_edge_labels = true;
			bool highlight_hair = true;
		};

		template <typename A, typename B, typename k, typename LabelFn>
		static std::string to_dot(
			const MetaGraph<A, B, k>& meta_graph,
			LabelFn&& label_fn,
			const Options& options
		) {
			std::ostringstream out;
			std::unordered_map<A, std::string> node_ids;
			std::unordered_set<A> hair_nodes(meta_graph.hair.begin(), meta_graph.hair.end());
			std::size_t next_id = 0;

			auto ensure_id = [&](const A& node) -> const std::string& {
				auto [it, inserted] = node_ids.try_emplace(node, "n" + std::to_string(next_id));
				if (inserted) {
					++next_id;
				}
				return it->second;
			};

			for (const auto& edge : meta_graph.edges) {
				ensure_id(edge.x);
				ensure_id(edge.y);
			}
			for (const auto& node : meta_graph.hair) {
				ensure_id(node);
			}

			out << "graph " << options.graph_name << " {\n";
			out << "  layout=neato;\n";
			out << "  overlap=false;\n";
			out << "  splines=true;\n";
			out << "  node [shape=circle, style=filled, fillcolor=white, fontname=\"Helvetica\"];\n";
			out << "  edge [fontname=\"Helvetica\"];\n";

			for (const auto& [node, node_id] : node_ids) {
				out << "  " << node_id << " [label=\"" << label_fn(node) << "\"";
				if (options.highlight_hair && hair_nodes.contains(node)) {
					out << ", fillcolor=\"#ffe082\"";
				}
				out << "];\n";
			}

			for (const auto& edge : meta_graph.edges) {
				out << "  " << node_ids.at(edge.x) << " -- " << node_ids.at(edge.y);
				if (options.show_edge_labels) {
					out << " [label=\"" << edge.common_boundary_components.size() << "\"]";
				}
				out << ";\n";
			}

			out << "}\n";
			return out.str();
		}

		template <typename A, typename B, typename k, typename LabelFn>
		static void write_dot(
			const std::filesystem::path& output_path,
			const MetaGraph<A, B, k>& meta_graph,
			LabelFn&& label_fn,
			const Options& options
		) {
			std::ofstream out(output_path);
			out << to_dot(meta_graph, std::forward<LabelFn>(label_fn), options);
		}
};
