"""Published amino-acid property tables (pure data, no third-party imports)."""

# Kyte & Doolittle (1982) hydropathy; Zamyatnin (1972) residue volume (A^3).
HYDROPATHY = dict(zip("ARNDCQEGHILKMFPSTWYV", [
    1.8, -4.5, -3.5, -3.5, 2.5, -3.5, -3.5, -0.4, -3.2, 4.5,
    3.8, -3.9, 1.9, 2.8, -1.6, -0.8, -0.7, -0.9, -1.3, 4.2], strict=True))
VOLUME = dict(zip("ARNDCQEGHILKMFPSTWYV", [
    88.6, 173.4, 114.1, 111.1, 108.5, 143.8, 138.4, 60.1, 153.2, 166.7,
    166.7, 168.6, 162.9, 189.9, 112.7, 89.0, 116.1, 227.8, 193.6, 140.0], strict=True))
# Tien et al. (2013) theoretical maximum accessible surface area (A^2).
MAX_ASA = dict(zip("ARNDCQEGHILKMFPSTWYV", [
    129.0, 274.0, 195.0, 193.0, 167.0, 225.0, 223.0, 104.0, 224.0, 197.0,
    201.0, 236.0, 224.0, 240.0, 159.0, 155.0, 172.0, 285.0, 263.0, 174.0], strict=True))
CHARGE_PH7 = {"D": -1, "E": -1, "K": 1, "R": 1}  # His treated as neutral (approximation)
