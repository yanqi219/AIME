#################################################
# Comparative Toxicogenomics Database AutoSearch
# Author: Qi Yan
# Date: 01/15/2019
#################################################
from bs4 import BeautifulSoup
import pandas as pd
from selenium import webdriver

########################
# PM 2.5
########################

# Set selenium environment for java
# download_dir = "C:/Users/QiYan/Dropbox/AIME/NCBI Search/CARB Scraping/PM2.5/"
options = webdriver.ChromeOptions()
options.add_experimental_option("prefs", {
  "download.default_directory": "C:/Users/QiYan/Dropbox/AIME/NCBI Search/CARB Scraping/PM2.5/",
  "download.prompt_for_download": False,
  "download.directory_upgrade": True,
  "safebrowsing.enabled": True
})
driver = webdriver.Chrome(options=options, executable_path='C:/Users/QiYan/Downloads/chromedriver_win32/chromedriver.exe')

# Begin scraping
year = ["2008","2009","2010","2011","2012","2013","2014","2015","2016"]
for i in range(0, len(year)):
    url = "http://www.arb.ca.gov/aqmis2/display.php?param=PM25&units=001&year="+year[i]+"&county_name=19-Los+Angeles&basin=--AIR+BASIN--&latitude=--PART+OF+STATE--&report=AQBYYR&order=basin%2Ccounty_name%2Cs.name&submit=Retrieve+Data&ptype=aqd&std15="

    # driver = webdriver.Chrome('C:/Users/QiYan/Downloads/chromedriver_win32/chromedriver.exe')
    # driver.implicitly_wait(30)
    driver.get(url)
    soup=BeautifulSoup(driver.page_source, 'lxml')

    hit_link = [a['href'] for a in soup.select('a')
                if a['href'].find('param=PM25') > -1]

    for j in range(0, len(hit_link)):
        station_url = 'https://www.arb.ca.gov/aqmis2/' + hit_link[j]
        driver.get(station_url)
        station_soup = BeautifulSoup(driver.page_source, 'lxml')   ### get monitor site info here!!!
        station_hit_link = [a['href'] for a in station_soup.select('a')
                            if a['href'].find('monitor=') > -1 & a['href'].find('DAVG') > -1]
        station_hit_link = list(set(station_hit_link))
        for k in range(0, len(station_hit_link)):
            monitor_url = 'https://www.arb.ca.gov/aqmis2/' + station_hit_link[k]
            driver.get(monitor_url)
            monitor_soup = BeautifulSoup(driver.page_source, 'lxml')
            monitor_hit_link = [a['href'] for a in monitor_soup.select('a')
                                if a['href'].find('download=y') > -1]
            download_url = 'https://www.arb.ca.gov/aqmis2/' + monitor_hit_link[0]   ### This is the URL!!! download it, get year, monitor site information
            driver.get(download_url)


########################
# Ozone
########################

# Set selenium environment for java

options = webdriver.ChromeOptions()
options.add_experimental_option("prefs", {
  "download.default_directory": "C:/Users/QiYan/Dropbox/AIME/NCBI Search/CARB Scraping/Ozone/",
  "download.prompt_for_download": False,
  "download.directory_upgrade": True,
  "safebrowsing.enabled": True
})
driver = webdriver.Chrome(options=options, executable_path='C:/Users/QiYan/Downloads/chromedriver_win32/chromedriver.exe')

# Begin scraping
year = ["2008","2009","2010","2011","2012","2013","2014","2015","2016"]
for i in range(0, len(year)):
    url = "https://www.arb.ca.gov/aqmis2/display.php?param=OZONE&units=007&year="+year[i]+"&county_name=19-Los+Angeles&basin=--AIR+BASIN--&latitude=--PART+OF+STATE--&report=AQBYYR&order=basin%2Ccounty_name%2Cs.name&submit=Retrieve+Data&ptype=aqd&std15="

    # driver = webdriver.Chrome('C:/Users/QiYan/Downloads/chromedriver_win32/chromedriver.exe')
    # driver.implicitly_wait(30)
    driver.get(url)
    soup=BeautifulSoup(driver.page_source, 'lxml')

    hit_link = [a['href'] for a in soup.select('a')
                if a['href'].find('param=OZONE') > -1]

    for j in range(0, len(hit_link)):
        monitor_url = 'https://www.arb.ca.gov/aqmis2/' + hit_link[j]
        driver.get(monitor_url)
        monitor_soup = BeautifulSoup(driver.page_source, 'lxml')
        monitor_hit_link = [a['href'] for a in monitor_soup.select('a')
                            if a['href'].find('download=y') > -1]
        download_url = 'https://www.arb.ca.gov/aqmis2/' + monitor_hit_link[0]  ### This is the URL!!! download it, get year, monitor site information
        driver.get(download_url)


########################
# NO2
########################

# Set selenium environment for java

options = webdriver.ChromeOptions()
options.add_experimental_option("prefs", {
  "download.default_directory": "C:/Users/QiYan/Dropbox/AIME/NCBI Search/CARB Scraping/NO2/",
  "download.prompt_for_download": False,
  "download.directory_upgrade": True,
  "safebrowsing.enabled": True
})
driver = webdriver.Chrome(options=options, executable_path='C:/Users/QiYan/Downloads/chromedriver_win32/chromedriver.exe')

# Begin scraping
year = ["2008","2009","2010","2011","2012","2013","2014","2015","2016"]
for i in range(0, len(year)):
    url = "https://www.arb.ca.gov/aqmis2/display.php?param=NO2&units=007&year="+year[i]+"&county_name=19-Los+Angeles&basin=--AIR+BASIN--&latitude=--PART+OF+STATE--&report=AQBYYR&order=basin%2Ccounty_name%2Cs.name&submit=Retrieve+Data&ptype=aqd&std15="

    # driver = webdriver.Chrome('C:/Users/QiYan/Downloads/chromedriver_win32/chromedriver.exe')
    # driver.implicitly_wait(30)
    driver.get(url)
    soup=BeautifulSoup(driver.page_source, 'lxml')

    hit_link = [a['href'] for a in soup.select('a')
                if a['href'].find('param=NO2') > -1]

    for j in range(0, len(hit_link)):
        monitor_url = 'https://www.arb.ca.gov/aqmis2/' + hit_link[j]
        driver.get(monitor_url)
        monitor_soup = BeautifulSoup(driver.page_source, 'lxml')
        monitor_hit_link = [a['href'] for a in monitor_soup.select('a')
                            if a['href'].find('download=y') > -1]
        download_url = 'https://www.arb.ca.gov/aqmis2/' + monitor_hit_link[0]  ### This is the URL!!! download it, get year, monitor site information
        driver.get(download_url)


########################
# NO
########################

# Set selenium environment for java

options = webdriver.ChromeOptions()
options.add_experimental_option("prefs", {
  "download.default_directory": "C:/Users/QiYan/Dropbox/AIME/NCBI Search/CARB Scraping/NO/",
  "download.prompt_for_download": False,
  "download.directory_upgrade": True,
  "safebrowsing.enabled": True
})
driver = webdriver.Chrome(options=options, executable_path='C:/Users/QiYan/Downloads/chromedriver_win32/chromedriver.exe')

# Begin scraping
year = ["2008","2009","2010","2011","2012","2013","2014","2015","2016"]
for i in range(0, len(year)):
    url = "https://www.arb.ca.gov/aqmis2/display.php?param=NO&units=007&year="+year[i]+"&county_name=19-Los+Angeles&basin=--AIR+BASIN--&latitude=--PART+OF+STATE--&report=AQBYYR&order=basin%2Ccounty_name%2Cs.name&submit=Retrieve+Data&ptype=aqd&std15="

    # driver = webdriver.Chrome('C:/Users/QiYan/Downloads/chromedriver_win32/chromedriver.exe')
    # driver.implicitly_wait(30)
    driver.get(url)
    soup=BeautifulSoup(driver.page_source, 'lxml')

    hit_link = [a['href'] for a in soup.select('a')
                if a['href'].find('param=NO&year') > -1]

    for j in range(0, len(hit_link)):
        monitor_url = 'https://www.arb.ca.gov/aqmis2/' + hit_link[j]
        driver.get(monitor_url)
        monitor_soup = BeautifulSoup(driver.page_source, 'lxml')
        monitor_hit_link = [a['href'] for a in monitor_soup.select('a')
                            if a['href'].find('download=y') > -1]
        download_url = 'https://www.arb.ca.gov/aqmis2/' + monitor_hit_link[0]  ### This is the URL!!! download it, get year, monitor site information
        driver.get(download_url)