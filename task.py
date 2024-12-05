import os
import math
import numpy as np
import geopandas as gpd
import psycopg2
import requests
import zipfile
from shapely.geometry import Point, Polygon, box
import folium
from pyproj import datadir
from shapely.geometry import shape

datadir.set_data_dir("/Users/Platon/anaconda3/share/proj")

# Константи
DATA_DIR = "data"
OVERPASS_URL = "https://overpass-api.de/api/interpreter"
SHAPEFILE_PATH = os.path.join(DATA_DIR, "ukraine_border.shp")
CONN_STRING = "dbname='test_Vodafone' user='karyna_b' password='1KARi061192' host='localhost' port='5433'"
RADIUS = 5 / 111  # Радіус 5 км у градусах
AZIMUTHS = [0, 120, 240]
ANGLE_WIDTH = 60


def get_ukraine_border(shapefile_path):
    """Отримує кордон України з локального файлу або Overpass API."""
    if os.path.exists(shapefile_path):
        # Завантаження з локального файлу
        ukraine = gpd.read_file(shapefile_path)
    else:
        # Запит до Overpass API
        query = """
        [out:json][timeout:25];
        (
          relation["name"="Україна"];
        );
        out body geom;
        """

        response = requests.get(OVERPASS_URL, params={"data": query})
        data = response.json()

        if "elements" not in data:
            raise ValueError("No data received from Overpass API!")

        for element in data["elements"]:
            if element["type"] == "relation" and "members" in element:
                coords = [
                    (member["geometry"][0]["lon"], member["geometry"][0]["lat"])
                    for member in element["members"]
                    if member["type"] == "way" and "geometry" in member
                ]
                if coords:
                    geometry = shape({"type": "Polygon", "coordinates": [coords]})
                    break
        else:
            raise ValueError("No border geometry found!")

        # Конвертація в GeoDataFrame
        ukraine = gpd.GeoDataFrame(
            [{"name": "Ukraine", "geometry": geometry}],
            crs="EPSG:4326"
        )

        # Збереження у форматі .shp
        os.makedirs(os.path.dirname(shapefile_path), exist_ok=True)
        ukraine.to_file(shapefile_path)

    return ukraine


def save_to_postgresql(geometry, table_name, conn_string):
    """Запис геометрії у PostgreSQL"""
    conn = psycopg2.connect(conn_string)
    cur = conn.cursor()
    cur.execute(f"""
        CREATE TABLE IF NOT EXISTS {table_name} (
            id SERIAL PRIMARY KEY,
            geom GEOMETRY(POLYGON, 4326)
        );
    """)
    for polygon in geometry:
        cur.execute(f"INSERT INTO {table_name} (geom) VALUES (ST_GeomFromText(%s, 4326));", (polygon.wkt,))
    conn.commit()
    cur.close()
    conn.close()

def create_grid(polygon, step_x=0.4, step_y=0.28):
    """Створення сітки квадратів з різними кроками по X та Y"""
    minx, miny, maxx, maxy = polygon.bounds
    x_coords = np.arange(minx, maxx, step_x)
    y_coords = np.arange(miny, maxy, step_y)
    grid = []
    for x in x_coords:
        for y in y_coords:
            cell = box(x, y, x + step_x, y + step_y)
            if polygon.intersects(cell):
                grid.append(cell)
    return grid

def save_grid_to_db(grid, conn_string, table_name="ukraine_grid"):
    """Збереження сітки у PostgreSQL"""
    conn = psycopg2.connect(conn_string)
    cur = conn.cursor()
    cur.execute(f"""
        CREATE TABLE IF NOT EXISTS {table_name} (
            id SERIAL PRIMARY KEY,
            geom GEOMETRY(POLYGON, 4326)
        );
    """)
    for cell in grid:
        cur.execute(f"INSERT INTO {table_name} (geom) VALUES (ST_GeomFromText(%s, 4326));", (cell.wkt,))
    conn.commit()
    cur.close()
    conn.close()

def create_sector(center, azimuth, radius, angle_width):
    """Створення сектора"""
    points = [center]
    start_angle = math.radians(azimuth - angle_width / 2)
    end_angle = math.radians(azimuth + angle_width / 2)
    for angle in np.arange(start_angle, end_angle, math.radians(1)):
        x = center.x + radius * math.cos(angle)
        y = center.y + radius * math.sin(angle)
        points.append(Point(x, y))
    points.append(center)
    return Polygon([[p.x, p.y] for p in points])

def save_intersection_to_db(intersections, conn_string, table_name="sector_intersections"):
    conn = psycopg2.connect(conn_string)
    cur = conn.cursor()
    cur.execute(f"""
        CREATE TABLE IF NOT EXISTS {table_name} (
            id SERIAL PRIMARY KEY,
            sector_id INT,
            square_id INT,
            intersects BOOLEAN
        );
    """)
    for row in intersections:
        cur.execute(f"""
            INSERT INTO {table_name} (sector_id, square_id, intersects)
            VALUES (%s, %s, %s);
        """, (row['sector_id'], row['square_id'], row['intersects']))
    conn.commit()
    cur.close()
    conn.close()
    
if __name__ == '__main__':
    # Завантаження кордону України
    ukraine = get_ukraine_border(SHAPEFILE_PATH)
    save_to_postgresql(ukraine.geometry.values, "ukraine_border", CONN_STRING)

    # Створення сітки
    ukraine_polygon = ukraine.geometry.unary_union
    grid = create_grid(ukraine_polygon,step_x=0.4, step_y=0.28)
    grid_gdf = gpd.GeoDataFrame({'geometry': grid}, crs=ukraine.crs)

    # Збереження сітки у файли та базу
    grid_gdf.to_file("ukraine_grid.geojson", driver="GeoJSON")
    save_grid_to_db(grid, CONN_STRING)

    # Створення секторів
    vertices = [Point(coord) for poly in grid for coord in poly.exterior.coords]
    sectors = [
        {"geometry": create_sector(vertex, azimuth, RADIUS, ANGLE_WIDTH)}
        for vertex in vertices for azimuth in AZIMUTHS
    ]
    sectors_gdf = gpd.GeoDataFrame(sectors, crs="EPSG:4326")

    # Обчислення перетинів
    intersections = []
    for square_id, square in enumerate(grid_gdf.geometry):
        for sector_id, sector in enumerate(sectors_gdf.geometry):
            if square.intersects(sector):
                intersections.append({
                    "sector_id": sector_id,
                    "square_id": square_id,
                    "intersects": True
                })
    save_intersection_to_db(intersections, CONN_STRING)

    # Візуалізація кордону, сітки та секторів
    map1 = folium.Map(location=[48.3794, 31.1656], zoom_start=6)
    folium.GeoJson(ukraine).add_to(map1)
    map1.save("ukraine_border.html")
    folium.GeoJson(grid_gdf).add_to(map1)
    map1.save("ukraine_with_grid.html")

    m = folium.Map(location=[48.3794, 31.1656], zoom_start=6)
    for _, row in grid_gdf.iterrows():
        folium.GeoJson(row['geometry'], style_function=lambda x: {'color': 'blue'}).add_to(m)
    for _, row in sectors_gdf.iterrows():
        folium.GeoJson(row['geometry'], style_function=lambda x: {'color': 'red', 'fillOpacity': 0.4}).add_to(m)
    m.save("ukraine_with_grid_and_sectors.html")