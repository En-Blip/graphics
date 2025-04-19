from PIL import Image
import csv

def create_image(red_file, green_file, blue_file):
    with open(red_file, 'r') as red, open(green_file, 'r') as green, open(blue_file, 'r') as blue:
        red_reader = csv.reader(red)
        green_reader = csv.reader(green)
        blue_reader = csv.reader(blue)

        # Read data from CSV files
        red_data = [list(map(int, row)) for row in red_reader]
        green_data = [list(map(int, row)) for row in green_reader]
        blue_data = [list(map(int, row)) for row in blue_reader]

    # Check if the dimensions of the CSV files match
    if len(red_data) != len(green_data) or len(green_data) != len(blue_data):
        print("Error: Dimensions of CSV files do not match.")
        sys.exit(1)

    # Create an image from the RGB data
    image = Image.new('RGB', (len(red_data[0]), len(red_data)))

    for y in range(len(red_data)):
        for x in range(len(red_data[0])):
            pixel = (red_data[y][x], green_data[y][x], blue_data[y][x])
            image.putpixel((x, y), pixel)

    # Save the image
    image.save('output_image.png')

if __name__ == "__main__":

    red_filename = 'imgR.csv';
    green_filename = 'imgG.csv';
    blue_filename = 'imgB.csv';
    create_image(red_filename, green_filename, blue_filename)
