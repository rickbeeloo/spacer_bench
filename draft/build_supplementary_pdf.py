import os
from md2pdf.core import md2pdf
import glob
from pypdf import PdfWriter, PdfReader



# Input and output paths
md_file = "supp/rawmd.notmd"
pdf_file = "../public/supplementary/supplementary.pdf"
style_file = "./supp/pdf_styles.css"
figures_dir = "supp/figures/supp/"  # Directory containing supplementary figure PDFs


md2pdf(pdf_file_path=pdf_file,
       md_file_path=md_file,
       css_file_path=style_file,
       base_url=os.path.dirname(md_file)
)

# # Read markdown content

# with open(md_file, 'r', encoding='utf-8') as f:
#     md_content = f.read()

# # Read CSS content
# with open(style_file, 'r', encoding='utf-8') as f:
#     css_content = f.read()

# # Split YAML header from content
# parts = md_content.split('---', 2)
# if len(parts) >= 3:
#     yaml_header = parts[1]
#     content = parts[2]
# else:
#     content = md_content

# # Process content to ensure proper heading hierarchy
# # Replace any ### that doesn't have a parent ## with ##
# lines = content.split('\n')
# processed_lines = []
# current_level = 0

# for line in lines:
#     if line.strip().startswith('#'):
#         level = len(re.match(r'^#+', line).group())
#         if level > current_level + 1:
#             # Reduce the heading level to maintain hierarchy
#             line = '#' * (current_level + 1) + line[level:]
#         current_level = len(re.match(r'^#+', line).group())
#     processed_lines.append(line)

# processed_content = '\n'.join(processed_lines)

# # Add the main content as a section
# pdf.add_section(
#     Section(
#         processed_content,
#         root=os.path.dirname(md_file),  # For correct image paths
#         paper_size="A4"
#     ),
#     user_css=css_content
# )

# # Set PDF metadata
# pdf.meta["title"] = title
# pdf.meta["creator"] = "markdown-pdf"

# # Save the PDF
# pdf.save(pdf_file)

print(f"Main supplementary PDF file created: {pdf_file}")

# def append_pdfs(base_pdf_path: str, pdfs_to_append: list[str], output_path: str) -> None:
#     """
#     Appends multiple PDFs to a base PDF file.
    
#     Args:
#         base_pdf_path: Path to the main PDF file
#         pdfs_to_append: List of paths to PDFs that should be appended
#         output_path: Path where the final merged PDF will be saved
#     """
#     pdf_merger = PdfWriter()

#     # Add the base PDF
#     with open(base_pdf_path, 'rb') as base_file:
#         pdf_reader = PdfReader(base_file)
#         pdf_merger.append(pdf_reader)

#     # Append each supplementary PDF
#     for pdf_path in pdfs_to_append:
#         print(f"Appending {pdf_path}")
#         with open(pdf_path, 'rb') as append_file:
#             pdf_reader = PdfReader(append_file)
#             pdf_merger.append(pdf_reader)

#     # Write the merged PDF
#     with open(output_path, 'wb') as output_file:
#         pdf_merger.write(output_file)

# # Find all supplementary figure PDFs
# supp_figure_pdfs = sorted(glob.glob(os.path.join(figures_dir, "*.pdf")))

# if supp_figure_pdfs:
#     # Create a temporary file for the merged PDF
#     temp_output = pdf_file + ".tmp"
    
#     # Merge PDFs
#     append_pdfs(pdf_file, supp_figure_pdfs, temp_output)
    
#     # Replace the original with the merged version
#     os.replace(temp_output, f"{pdf_file.replace('.pdf', '')}_with_supplements.pdf")
#     print(f"Appended {len(supp_figure_pdfs)} supplementary figure PDFs to {pdf_file}")
#     print("Appended PDFs:", "\n".join(supp_figure_pdfs))
#     print(f"Main supplementary PDF file created: {pdf_file.replace('.pdf', '')}_with_supplements.pdf")
# else:
#     print("No supplementary figure PDFs found in", figures_dir)