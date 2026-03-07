//! Available E2 measurement units.
//!
//! Mirrors the `availableUnits` object in modeldata.js.

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Unit {
    PgPerMl,
    PmolPerL,
    NgPerL,
    /// Easter egg from the original: firkin per furlong³
    Fff,
}

#[derive(Debug, Clone, Copy)]
pub struct UnitInfo {
    pub label: &'static str,
    pub conversion_factor: f64,
    /// Decimal places to show in output
    pub precision: usize,
}

impl Unit {
    pub fn info(self) -> UnitInfo {
        match self {
            Self::PgPerMl => UnitInfo {
                label: "pg/mL",
                conversion_factor: 1.0,
                precision: 0,
            },
            Self::PmolPerL => UnitInfo {
                label: "pmol/L",
                conversion_factor: 3.6713,
                precision: 0,
            },
            Self::NgPerL => UnitInfo {
                label: "ng/L",
                conversion_factor: 1.0,
                precision: 0,
            },
            Self::Fff => UnitInfo {
                label: "firkin/furlong\u{B3}",
                conversion_factor: 0.000320496,
                precision: 4,
            },
        }
    }

    pub fn from_str(s: &str) -> Option<Self> {
        match s {
            "pg/mL" => Some(Self::PgPerMl),
            "pmol/L" => Some(Self::PmolPerL),
            "ng/L" => Some(Self::NgPerL),
            "FFF" => Some(Self::Fff),
            _ => None,
        }
    }

    pub fn conversion_factor(self) -> f64 {
        self.info().conversion_factor
    }
}
