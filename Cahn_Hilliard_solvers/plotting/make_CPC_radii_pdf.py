from reportlab.pdfgen import canvas
from reportlab.lib.pagesizes import letter, landscape

from PIL import Image
import os
from glob import glob
import re

#unused
def add_images_to_pdf(image_paths, pdf_path):
    c = canvas.Canvas(pdf_path, pagesize=letter)
    
    for image_path in image_paths:
        # Open the image using PIL
        img = Image.open(image_path)
        
        # Get the image size
        img_width, img_height = img.size
        
        # Calculate the position and size of the image in the PDF
        # Scale the image to fit the page width, maintaining aspect ratio
        aspect = img_width / float(img_height)
        if aspect > 1:
            width = 550
            height = 550 / aspect
        else:
            height = 700
            width = 700 * aspect
        
        # Add the image to the PDF
        c.drawImage(image_path, 30, 750 - height, width=width, height=height)
        
        # Add a new page for the next image
        c.showPage()
    
    c.save()
#unused
def add_images_to_pdf_grid(image_paths, pdf_path, images_per_row=2):
    c = canvas.Canvas(pdf_path, pagesize=letter)
    width, height = letter

    # Dimensions for the grid
    margin = 30
    padding = 10
    img_width = (width - 2 * margin - (images_per_row - 1) * padding) / images_per_row

    current_x = margin
    current_y = height - margin

    for index, image_path in enumerate(image_paths):
        # Open the image using PIL
        img = Image.open(image_path)

        # Maintain aspect ratio
        aspect = img.width / float(img.height)
        img_height = img_width / aspect

        # Adjust the y position if the image is too tall
        if current_y - img_height < margin:
            c.showPage()
            current_y = height - margin

        # Draw the image
        c.drawImage(image_path, current_x, current_y - img_height, width=img_width, height=img_height)
        
        # Update the position for the next image
        current_x += img_width + padding
        if (index + 1) % images_per_row == 0:
            current_x = margin
            current_y -= img_height + padding

    c.save()

def extract_variables(image_path):
    # Extract CPC and cohesin numbers using regex
    cpc_match = re.search(r'CPC_(\d+)', image_path)
    cohesin_match = re.search(r'cohesin_(\d+)', image_path)
    cpc = int(cpc_match.group(1)) if cpc_match else 0
    cohesin = int(cohesin_match.group(1)) if cohesin_match else 0
    return cpc, cohesin
def draw_rotated_text(c, text, x, y, angle):
    c.saveState()
    c.translate(x, y)
    c.rotate(angle)
    c.drawString(0, 0, text)
    c.restoreState()
def add_images_to_pdf_sorted_grid(image_paths, pdf_path):
    # Extract CPC and cohesin numbers and sort based on them
    images_info = [(extract_variables(path), path) for path in image_paths]
    images_info.sort(key=lambda x: (x[0][1], x[0][0]))  # Sort by cohesin, then CPC

    # Create a set of unique CPC and cohesin values
    unique_cpc = sorted({extract_variables(path)[0] for path in image_paths})
    unique_cohesin = sorted({extract_variables(path)[1] for path in image_paths})

    c = canvas.Canvas(pdf_path, pagesize=landscape(letter))
    page_width, page_height = landscape(letter)

      # Create a set of unique CPC and cohesin values
    unique_cpc = sorted({extract_variables(path)[0] for path in image_paths})
    unique_cohesin = sorted({extract_variables(path)[1] for path in image_paths})


    # Dimensions for the grid
    margin = 30
    padding = 10
    img_width = (page_width - 2 * margin - (len(unique_cpc) - 1) * padding) / len(unique_cpc)

    # Calculate the positions for each unique CPC and cohesin value
    cpc_to_x = {cpc: margin + i * (img_width + padding) for i, cpc in enumerate(unique_cpc)}
    cohesin_to_y = {cohesin: page_height - margin - (i * (img_width + padding)) for i, cohesin in enumerate(unique_cohesin)}

    # Add labels for CPC columns
    c.drawString(page_width/2 -20, page_height - 10, f'CPC radius')

    for cpc, x in cpc_to_x.items():
        c.drawString(x + img_width / 2, page_height - 25, f'{cpc/256}')


    # Add labels for cohesin rows
    draw_rotated_text(c, "(Full) cohesin width", margin -15, (page_height - margin) / 2 - 15, 90)

    for cohesin, y in cohesin_to_y.items():
        # c.drawString(margin / 2, y - img_width / 2, f'cohesin {cohesin}')
        c.saveState()
        c.translate(25, y - img_width / 2 +10)
        c.rotate(90)
        c.drawString(0, 0, f'{(cohesin*2)/256}')
        c.restoreState()
        # Draw x-axis (CPC values)
    # c.setStrokeColorRGB(0, 0, 0)
    # c.setLineWidth(1)
    # c.line(margin, margin, page_width - margin, margin)
    # for cpc, x in cpc_to_x.items():
    #     c.drawString(x, margin - 20, str(cpc))

    # # Draw y-axis (cohesin values)
    # c.line(margin, margin, margin, page_height - margin)
    # for cohesin, y in cohesin_to_y.items():
    #     c.drawString(margin - 40, y, str(cohesin))


    for (cpc, cohesin), image_path in images_info:
        current_x = cpc_to_x[cpc]
        current_y = cohesin_to_y[cohesin]

        # Open the image using PIL
        img = Image.open(image_path)

        # Maintain aspect ratio
        aspect = img.width / float(img.height)
        img_height = img_width / aspect

        # Draw the image
        c.drawImage(image_path, current_x, current_y - img_height, width=img_width, height=img_height)

    c.save()

indir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/plotting/"
paths = glob(f"{indir}/CPC_radii/*/")
image_paths = []
epsilon = "0.08"
for p in paths:
    print(p)
    eps_match = re.search(r'eps_(\d*\.?\d*)', p)
    eps = (eps_match.group(1)) if eps_match else 0
    if eps == epsilon:
        png = f"{p}radii_over_time.png"
        image_paths.append(png)
print(image_paths)
pdf_path = f"{indir}/radii_over_time_{epsilon}.pdf"
# add_images_to_pdf(image_paths, pdf_path)

add_images_to_pdf_sorted_grid(image_paths, pdf_path)
