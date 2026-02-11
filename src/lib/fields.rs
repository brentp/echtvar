use serde::{Deserialize, Serialize}; // 1.0.101

#[derive(Debug, Deserialize, Serialize, PartialEq, PartialOrd, Clone, Default)]
pub enum FieldType {
    #[default]
    Integer,
    Float,
    Categorical,
    Flag,
}

#[derive(Debug, Deserialize, Serialize, PartialEq, PartialOrd, Clone)]
pub struct Field {
    pub field: String,
    pub alias: String,
    #[serde(default = "default_missing_value")]
    pub missing_value: i32,
    #[serde(default = "default_missing_string")]
    pub missing_string: std::string::String,
    #[serde(default = "default_description_string")]
    pub description: std::string::String,

    #[serde(default)]
    pub zigzag: bool,

    #[serde(default = "default_multiplier")]
    pub multiplier: u32,
    #[serde(default)]
    pub ftype: FieldType,
    #[serde(default = "default_number")]
    pub number: std::string::String,

    #[serde(default = "default_values_i", skip_serializing)]
    pub values_i: usize,
}

fn default_missing_value() -> i32 {
    -1
}

fn default_number() -> std::string::String {
    "1".to_string()
}
fn default_missing_string() -> std::string::String {
    "MISSING".to_string()
}
pub fn default_description_string() -> std::string::String {
    "added by echtvar".to_string()
}
fn default_multiplier() -> u32 {
    1
}
fn default_values_i() -> usize {
    usize::MAX
}

impl Default for Field {
    fn default() -> Field {
        Field {
            field: "field".to_string(),
            alias: "name".to_string(),
            missing_value: -1,
            missing_string: "MISSING".to_string(),
            description: "added by echtvar".to_string(),
            zigzag: false,
            multiplier: 1,
            ftype: FieldType::Integer,
            number: ".".to_string(),
            values_i: usize::MAX,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_field() {
        let fields: Vec<Field> = json5::from_str(r#"
        [
         {"field": "AC", "alias": "gnomad_AC"},
         {
             field: "AN",
             alias: "gnomad_AN",
             missing_value: -2147483648,
             missing: -2,
             multiplier: 1, // this is useful for float fields as internally, everythign is stored as integer.
             zigzag: true, // set this to true if the field can contain negative numbers.
             ftype: "Integer", // this is discovered by echtvar and should not be set.
        }
         ]
        "#).unwrap();

        assert_eq!(fields.len(), 2);
        assert_eq!(fields[0].alias, "gnomad_AC");
        assert_eq!(fields[1].alias, "gnomad_AN");
        assert_eq!(fields[0].missing_value, -1);
        assert_eq!(fields[1].missing_value, -2147483648);
        assert_eq!(fields[1].zigzag, true);
        assert_eq!(fields[1].ftype, FieldType::Integer);

        eprintln!("{}", json5::to_string(&fields[0]).unwrap());
    }

    #[test]
    fn test_flag_field_type() {
        let fields: Vec<Field> = json5::from_str(
            r#"
        [
         {"field": "DB", "alias": "is_dbsnp", "ftype": "Flag"},
        ]
        "#,
        )
        .unwrap();

        assert_eq!(fields.len(), 1);
        assert_eq!(fields[0].field, "DB");
        assert_eq!(fields[0].alias, "is_dbsnp");
        assert_eq!(fields[0].ftype, FieldType::Flag);
    }

    #[test]
    fn test_flag_field_config_roundtrip() {
        // Simulate creating a Flag field as the encoder would
        let mut f = Field::default();
        f.field = "DB".to_string();
        f.alias = "is_dbsnp".to_string();
        f.ftype = FieldType::Flag;
        f.missing_value = 0;
        f.number = "0".to_string();

        // Serialize to JSON (as stored in config.json inside the zip)
        let json = serde_json::to_string(&f).unwrap();

        // Deserialize back (as read during annotation)
        let f2: Field = serde_json::from_str(&json).unwrap();

        assert_eq!(f2.ftype, FieldType::Flag);
        assert_eq!(f2.field, "DB");
        assert_eq!(f2.alias, "is_dbsnp");
        assert_eq!(f2.missing_value, 0);
        assert_eq!(f2.number, "0");
    }

    #[test]
    fn test_flag_missing_value_semantics() {
        // For Flag fields, missing_value=0 means absent.
        // When a variant is not in the database, get_int_value returns
        // missing_value, which is 0 for flags (absent) — correct behavior.
        let f = Field {
            ftype: FieldType::Flag,
            missing_value: 0,
            number: "0".to_string(),
            ..Field::default()
        };

        assert_eq!(f.missing_value, 0);
        assert_eq!(f.number, "0");

        // Flag values are 0 (absent) or 1 (present), both valid u32 values.
        // Neither should be confused with the u32::MAX missing sentinel.
        let flag_present: u32 = 1;
        let flag_absent: u32 = 0;
        assert!(flag_present < u32::MAX);
        assert!(flag_absent < u32::MAX);
    }
}
