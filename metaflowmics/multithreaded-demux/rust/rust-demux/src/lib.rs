#[derive(Clone, Debug)]
pub struct Barcode {
	pub sample: String,
	pub fwd: Vec<u8>,
	pub rev: Vec<u8>,
}
