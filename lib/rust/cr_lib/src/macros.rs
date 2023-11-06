#[macro_export]
macro_rules! per_type_metric {
    ($new_metric_name:ident, $new_visitor_name:ident, $key_type:ident, $inner_type:ident) => {
        #[derive(Clone, Deserialize, Serialize, Metric)]
        pub struct $new_metric_name(pub metric::TxHashMap<$key_type, $inner_type>);

        pub struct $new_visitor_name {
            pub metrics: $new_metric_name,
        }

        impl $new_visitor_name {
            pub fn new() -> Self {
                $new_visitor_name {
                    metrics: metric::Metric::new(),
                }
            }
        }

        impl Default for $new_visitor_name {
            fn default() -> Self {
                $new_visitor_name::new()
            }
        }
    };
}
