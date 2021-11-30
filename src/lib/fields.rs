use serde::{Deserialize, Serialize}; // 1.0.101

#[derive(Debug, Deserialize, Serialize, PartialEq)]
pub enum FieldType {
    Integer,
    Float,
    Categorical,
}

#[derive(Debug, Deserialize, Serialize)]
pub struct Field {
    pub field: String,
    pub alias: String,
    #[serde(default = "default_missing_value")]
    pub missing_value: i32,
    #[serde(default)]
    pub zigzag: bool,
    #[serde(default = "default_multiplier")]
    pub multiplier: u32,
    #[serde(default)]
    pub ftype: FieldType,
}

fn default_missing_value() -> i32 {
    -1
}
fn default_multiplier() -> u32 {
    1
}

impl Default for Field {
    fn default() -> Field {
        Field {
            field: "field".to_string(),
            alias: "name".to_string(),
            missing_value: -1,
            zigzag: false,
            multiplier: 1,
            ftype: FieldType::Integer,
        }
    }
}

impl Default for FieldType {
    fn default() -> FieldType {
        FieldType::Integer
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_read() {
        let fields: Vec<Field> = json5::from_str(r#"
        [
         {"field": "AC", "alias": "gnomad_AC"},
         {
             field: "AN",
             alias: "gnomad_AN",
             missing_value: -2147483648,
             multiplier: 1, // this is useful for float fields as internally, everythign is stored as integer.
             zigzag: false, // set this to true if the field can contain negative numbers.
             ftype: "Integer", // this is discovered by echtvar and should not be set.
        }
         ]
        "#).unwrap();

        assert_eq!(fields.len(), 2);
        assert_eq!(fields[0].alias, "gnomad_AC");
        assert_eq!(fields[1].alias, "gnomad_AN");
        assert_eq!(fields[0].missing_value, -1);
        assert_eq!(fields[1].missing_value, -2147483648);
        assert_eq!(fields[1].ftype, FieldType::Integer);

        eprintln!("{}", json5::to_string(&fields[0]).unwrap());
    }
}
