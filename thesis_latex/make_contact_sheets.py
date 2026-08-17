from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


ROOT = Path(__file__).resolve().parent
PAGES = ROOT / "tmp" / "pdf_final" / "pages"
SHEETS = ROOT / "tmp" / "pdf_final" / "contact_sheets"


def main() -> None:
    SHEETS.mkdir(parents=True, exist_ok=True)
    pages = sorted(PAGES.glob("page-*.png"))
    if not pages:
        raise RuntimeError("No rendered PDF pages found")

    font = ImageFont.load_default()
    for offset in range(0, len(pages), 4):
        batch = pages[offset : offset + 4]
        opened = [Image.open(path).convert("RGB") for path in batch]
        thumb_width = 720
        thumbnails = []
        for page in opened:
            ratio = thumb_width / page.width
            thumbnails.append(
                page.resize(
                    (thumb_width, round(page.height * ratio)),
                    Image.Resampling.LANCZOS,
                )
            )
        cell_height = max(image.height for image in thumbnails) + 28
        sheet = Image.new("RGB", (thumb_width * 2, cell_height * 2), "white")
        draw = ImageDraw.Draw(sheet)
        for index, (image, source) in enumerate(zip(thumbnails, batch)):
            x = (index % 2) * thumb_width
            y = (index // 2) * cell_height
            page_number = int(source.stem.split("-")[-1])
            draw.text((8 + x, 6 + y), f"PDF page {page_number}", fill="black", font=font)
            sheet.paste(image, (x, y + 24))
        first = offset + 1
        last = offset + len(batch)
        sheet.save(SHEETS / f"pages-{first:03d}-{last:03d}.png")

    print(f"Pages: {len(pages)}")
    print(f"Contact sheets: {(len(pages) + 3) // 4}")


if __name__ == "__main__":
    main()
