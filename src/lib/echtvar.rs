use crate::var32;

pub struct EchtVar<T> {
    pub name: String,
    pub missing: T,
    pub values: Vec<T>,
}

pub struct EchtVars<'a> {
    pub zip: &'a mut zip::ZipArchive<std::fs::File>,
    pub chrom: String,
    pub start: u32,
    pub var32s: Vec<var32::Var32>,
    pub longs: Vec<var32::LongVariant>,
    pub ints: Vec<EchtVar<i32>>,
    pub floats: Vec<EchtVar<f32>>,
}


