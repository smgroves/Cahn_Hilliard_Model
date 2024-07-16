from reportlab.pdfgen import canvas
from reportlab.lib.pagesizes import letter, landscape

from PIL import Image
import os
from glob import glob
import re

def extract_variables(image_path):
    # Extract CPC and alpha numbers using regex
    eps_match = re.search(r'eps_(\d*\.?\d*)', image_path) #updated for um
    alpha_match = re.search(r'alpha_(-?\d*\.?\d*)', image_path) #updated for um

    eps = float(eps_match.group(1)) if eps_match else 0
    alpha = float(alpha_match.group(1)) if alpha_match else 0
    return eps, alpha
def draw_rotated_text(c, text, x, y, angle):
    c.saveState()
    c.translate(x, y)
    c.rotate(angle)
    c.drawString(0, 0, text)
    c.restoreState()
def add_images_to_pdf_sorted_grid(image_paths, pdf_path, title):
    # Extract CPC and alpha numbers and sort based on them
    images_info = [(extract_variables(path), path) for path in image_paths]
    images_info.sort(key=lambda x: (x[0][1], x[0][0]))  # Sort by eps, then alpha

    # Create a set of unique CPC and alpha values
    unique_eps = sorted({extract_variables(path)[0] for path in image_paths})
    unique_alpha = sorted({extract_variables(path)[1] for path in image_paths})

    c = canvas.Canvas(pdf_path, pagesize=landscape(letter))
    page_width, page_height = landscape(letter)
    # page_width = 792
    # page_height = 792
    # print(page_height, page_width)
    # Dimensions for the grid
    margin = 30
    padding = 0
    img_width = (page_width - 2 * margin - (len(unique_eps) - 1) * padding) / len(unique_eps)

    # Calculate the positions for each unique CPC and alpha value
    eps_to_x = {eps: margin + i * (img_width + padding) for i, eps in enumerate(unique_eps)}
    alpha_to_y = {alpha: page_height - margin - (i * (img_width + padding)) for i, alpha in enumerate(unique_alpha)}
    print(alpha_to_y)
    c.drawString(10, page_height - 20, title)

    # Add labels for CPC columns
    c.drawString(page_width/2 -10, 7, f'Epsilon')

    for eps, x in eps_to_x.items():
        c.drawString(x + img_width / 2 - 20,  20, f'{eps}')


    # Add labels for alpha rows
    draw_rotated_text(c, "Alpha", margin -15, (page_height - margin) / 2 - 15, 90)

    for alpha, y in alpha_to_y.items():
        # c.drawString(margin / 2, y - img_width / 2, f'alpha {alpha}')
        c.saveState()
        c.translate(25, y - img_width / 2 +10)
        c.rotate(90)
        c.drawString(0, 0, f'{(alpha)}')
        c.restoreState()
        # Draw x-axis (CPC values)
    # c.setStrokeColorRGB(0, 0, 0)
    # c.setLineWidth(1)
    # c.line(margin, margin, page_width - margin, margin)
    # for eps, x in eps_to_x.items():
    #     c.drawString(x, margin - 20, str(eps))

    # # Draw y-axis (alpha values)
    # c.line(margin, margin, margin, page_height - margin)
    # for alpha, y in alpha_to_y.items():
    #     c.drawString(margin - 40, y, str(alpha))


    for (eps, alpha), image_path in images_info:
        current_x = eps_to_x[eps]
        current_y = alpha_to_y[alpha]

        # Open the image using PIL
        img = Image.open(image_path)

        # Maintain aspect ratio
        aspect = img.width / float(img.height)
        img_height = img_width / aspect

        # Draw the image
        c.drawImage(image_path, current_x, current_y - img_height, width=img_width, height=img_height)

    c.save()

time = "0"
# sim=f"03_25_24_CPC_relaxed_RefModel_128x64_03_25_24_relaxed_RefModel_{time}_256x256_{time}s_8.4max"
# sim=f"10_24_23_CPC_tensed_RefModel_128x64_post_transition_07_14_24_500s_post_transition_base_20Pac_{time}_256x256_{time}s_8.4max"
sim=f"10_16_23_CPC_relaxed_RefModel_128x64_10_16_23_relaxed_RefModel_Mps1_phos_Plk1a_20Pac_transactiv_70_256x256_70s_8.4max"
# sim = f"10_16_23_CPC_relaxed_RefModel_128x64_10_16_23_relaxed_RefModel_Mps1_phos_Plk1a_20Pac_transactiv_{time}_256x256_{time}s_8.4max"
indir = f"/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/plotting/radii_over_time_level_set_plots/{sim}/"
paths = glob(f"{indir}/*/")
image_paths = []
for p in paths:
    png = f"{p}kymograph_x_128.png"
    image_paths.append(png)

print(image_paths)
pdf_path = f"/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/plotting/kymograph_pdfs/{sim}.pdf"
# add_images_to_pdf(image_paths, pdf_path)

add_images_to_pdf_sorted_grid(image_paths, pdf_path, title = sim)
