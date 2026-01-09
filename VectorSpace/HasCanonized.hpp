#include <type_traits>
#include <concepts>

namespace VectorSpace {

template<typename T, typename k>
concept HasCanonized =
    std::is_class_v<T> &&
    requires (BasisElement<T,k>& e) {
        { T::canonized(e) } -> std::same_as<BasisElement<T,k>>;
    };

} //
