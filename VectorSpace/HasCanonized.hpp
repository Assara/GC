template<typename T, typename k>
concept HasCanonized =
    requires (BasisElement<T,k> e) {
        { T::canonized(e) } -> std::same_as<BasisElement<T,k>>;
    };
