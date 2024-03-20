#[macro_export]
macro_rules! per_type_metric {
    ($new_metric_name:ident, $new_visitor_name:ident, $key_type:ident, $inner_type:ident) => {
        #[derive(Default, Serialize, Deserialize, Metric)]
        pub struct $new_metric_name(pub metric::TxHashMap<$key_type, $inner_type>);

        pub struct $new_visitor_name {
            pub metrics: $new_metric_name,
        }

        impl $new_visitor_name {
            pub fn new() -> Self {
                $new_visitor_name {
                    metrics: Default::default(),
                }
            }
        }
    };
}
