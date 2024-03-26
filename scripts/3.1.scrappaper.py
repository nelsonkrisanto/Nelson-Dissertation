'''
Forked from 'ScrapPaper' by M. R. Rafsanjani.
Updated script to merge scraped Pubmed Summary / Abstract format data.
e.g.: https://pubmed.ncbi.nlm.nih.gov/?term=dengue+virus+primer
'''

import requests
import random
import time
import pandas as pd
from bs4 import BeautifulSoup
import os
from fuzzywuzzy import process, fuzz
import re

# Initialize
print("Initiating... please wait.\n")

# Define functions
def wait():
    time.sleep(random.uniform(1, 3))  # Use a smaller range to make it slightly faster
    print("Continuing...\n")

def get_pubmed_data(URL, headers, is_abstract_format=False):
    data = []
    page_num, page_view = 1, 200
    page_total_num = 1  # Default to 1 in case the first request fails

    try:
        URL_edit = URL + "&page=" + str(page_num) + "&size=" + str(page_view)
        page = requests.get(URL_edit, headers=headers)
        soup = BeautifulSoup(page.content, "html.parser")
        wait()

        # Get total number of pages
        page_total = soup.find("label", class_="of-total-pages").text
        page_total_num = int(''.join(filter(str.isdigit, page_total)))

    except Exception as e:
        print(f"Error: {e}")
        return data

    for i in range(page_total_num):
        URL_edit = URL + "&page=" + str(i + 1) + "&size=" + str(page_view)
        page = requests.get(URL_edit, headers=headers)
        soup = BeautifulSoup(page.content, "html.parser")
        wait()

        if is_abstract_format:
            results = soup.find_all("div", class_="results-article")
            for result in results:
                title = result.find("h1", class_="heading-title").get_text(strip=True)
                reference = result.find("span", class_="cit").get_text(strip=True) if result.find("span", class_="cit") else ""
                pubmed_link = "https://pubmed.ncbi.nlm.nih.gov/" + result.find("strong", class_="current-id").get_text(strip=True)
                full_text_link = result.find("div", class_="full-text-links-list").find("a", href=True)['href'] if result.find("div", class_="full-text-links-list") else ""
                data.append([title, reference, pubmed_link, full_text_link])
        else:
            results = soup.find_all("article", class_="full-docsum")
            for result in results:
                title = result.find("a", class_="docsum-title").get_text(strip=True)
                reference = result.find("span", class_="docsum-journal-citation full-journal-citation").get_text(strip=True)
                pubmed_link = "https://pubmed.ncbi.nlm.nih.gov" + result.find("a")['href']
                data.append([title, reference, pubmed_link, ""])

    return data

def merge_data(data1, data2):
    merged_data = []
    for row1 in data1:
        title1, reference1, pubmed_link1, _ = row1
        match = process.extractOne(title1, [row2[0] for row2 in data2], scorer=fuzz.partial_ratio)
        if match and match[1] >= 80:  # If similarity is 80% or more
            matched_row = data2[[row2[0] for row2 in data2].index(match[0])]
            _, _, _, full_text_link2 = matched_row
            merged_data.append([title1, reference1, pubmed_link1, full_text_link2])
        else:
            merged_data.append(row1)
    return merged_data

# Define a function to extract the year
def extract_year(reference):
    # Regex pattern to find the year
    year_pattern = re.compile(r'\b(\d{4})\b')
    year_match = year_pattern.search(reference)
    if year_match:
        return year_match.group(1)
    else:
        return "Year not found"

# Getting and setting the URL
URL_input = input("Please paste search URL and press Enter: ")
headers = {
    'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/119.0.0.0 Safari/537.36'
}

# Extract data
data_summary = get_pubmed_data(URL_input, headers)
data_abstract = get_pubmed_data(URL_input, headers, is_abstract_format=True)

# Merge data based on title similarity
merged_data = merge_data(data_summary, data_abstract)

# Modify the final data extraction to include year
final_data = []
for title, ref, pubmed, full_text in merged_data:
    year = extract_year(ref)
    final_data.append([title, year, ref, pubmed, full_text])

# Convert to DataFrame and save to Excel
df = pd.DataFrame(final_data, columns=['Title', 'Year', 'Reference', 'Pubmed Link', 'Full Text Link'])
df.to_excel('merged_pubmed_data.xlsx', index=False)

print("Job finished, data merged and saved to 'merged_pubmed_data.xlsx'.")
