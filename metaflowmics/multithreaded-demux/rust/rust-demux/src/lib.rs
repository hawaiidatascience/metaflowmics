// #[derive(Copy, Clone, Debug, PartialEq)]
// pub enum Dna {
// 	A = 65,
// 	C = 67,
// 	G = 71,
// 	T = 84,
// 	N = 78
// }

// impl From<u8> for Dna {
//     fn from(byte: u8) -> Self {
// 		match byte {
// 			65 => Dna::A,
// 			67 => Dna::C,
// 			71 => Dna::G,
// 			84 => Dna::T,
// 			78 => Dna::N,
// 			_ => panic!("Invalid nucleotide")
// 		}
//     }
// }

#[derive(Clone, Debug)]
pub struct Barcode {
	pub sample: String,
	pub fwd: Vec<u8>,
	pub rev: Vec<u8>,
}
