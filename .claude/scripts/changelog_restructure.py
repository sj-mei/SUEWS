#!/usr/bin/env python3
"""
Comprehensive CHANGELOG restructuring script that:
1. Parses CHANGELOG.md to structured JSON
2. Cleans duplicates and merges entries
3. Sorts properly by date
4. Formats back to clean markdown
"""

import re
import json
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Tuple
import hashlib


class ChangelogRestructurer:
    def __init__(self, changelog_path: str = 'CHANGELOG.md'):
        self.changelog_path = Path(changelog_path)
        self.date_pattern = re.compile(r'^- (\d{1,2} \w+ \d{4}):\s*$')
        self.entry_pattern = re.compile(r'^  - \[(\w+)\] (.+)$')
        
    def parse_to_json(self) -> Dict:
        """Parse CHANGELOG.md into structured JSON format"""
        if not self.changelog_path.exists():
            raise FileNotFoundError(f"{self.changelog_path} not found")
            
        content = self.changelog_path.read_text()
        lines = content.split('\n')
        
        # Structure to hold parsed data
        data = {
            'header': [],
            'entries': {},  # date_str -> list of entries
            'raw_entries': []  # For debugging
        }
        
        current_date = None
        current_entries = []
        in_header = True
        
        i = 0
        while i < len(lines):
            line = lines[i]
            
            # Check for date line
            date_match = self.date_pattern.match(line)
            if date_match:
                # Save previous date's entries
                if current_date and current_entries:
                    if current_date not in data['entries']:
                        data['entries'][current_date] = []
                    data['entries'][current_date].extend(current_entries)
                    data['raw_entries'].append({
                        'date': current_date,
                        'entries': current_entries.copy()
                    })
                
                # Start new date
                current_date = date_match.group(1)
                current_entries = []
                in_header = False
                
            elif in_header:
                data['header'].append(line)
                
            elif line.startswith('  - '):
                # This is an entry line
                entry_match = self.entry_pattern.match(line)
                if entry_match:
                    category = entry_match.group(1)
                    description = entry_match.group(2)
                    
                    # Start new entry
                    current_entry = {
                        'category': category,
                        'description': description,
                        'details': []
                    }
                    
                    # Look ahead for detail lines (starting with '    -')
                    j = i + 1
                    while j < len(lines) and lines[j].startswith('    -'):
                        detail = lines[j].strip()[2:]  # Remove '- '
                        current_entry['details'].append(detail)
                        j += 1
                    
                    current_entries.append(current_entry)
                    i = j - 1  # Skip the detail lines we just processed
                    
            i += 1
        
        # Don't forget the last entry
        if current_date and current_entries:
            if current_date not in data['entries']:
                data['entries'][current_date] = []
            data['entries'][current_date].extend(current_entries)
            data['raw_entries'].append({
                'date': current_date,
                'entries': current_entries
            })
        
        return data
    
    def clean_entries(self, data: Dict) -> Dict:
        """Clean duplicates and merge similar entries"""
        cleaned_data = data.copy()
        cleaned_entries = {}
        
        # Track seen entries to remove duplicates
        seen_entries = set()
        
        for date_str, entries in data['entries'].items():
            cleaned_date_entries = []
            
            for entry in entries:
                # Create a hash of the entry to detect duplicates
                entry_hash = hashlib.md5(
                    json.dumps(entry, sort_keys=True).encode()
                ).hexdigest()
                
                if entry_hash not in seen_entries:
                    seen_entries.add(entry_hash)
                    cleaned_date_entries.append(entry)
            
            if cleaned_date_entries:
                cleaned_entries[date_str] = cleaned_date_entries
        
        cleaned_data['entries'] = cleaned_entries
        return cleaned_data
    
    def sort_entries(self, data: Dict) -> List[Tuple[str, List]]:
        """Sort entries by date in reverse chronological order"""
        date_entries = []
        
        for date_str, entries in data['entries'].items():
            try:
                date_obj = datetime.strptime(date_str, '%d %b %Y')
                date_entries.append((date_obj, date_str, entries))
            except ValueError:
                print(f"Warning: Could not parse date '{date_str}'")
        
        # Sort by date object in reverse order
        date_entries.sort(key=lambda x: x[0], reverse=True)
        
        # Return sorted list of (date_str, entries) tuples
        return [(date_str, entries) for _, date_str, entries in date_entries]
    
    def format_to_markdown(self, data: Dict, sorted_entries: List[Tuple[str, List]]) -> str:
        """Format the cleaned and sorted data back to markdown"""
        lines = []
        
        # Add header
        lines.extend(data['header'])
        
        # Add sorted entries
        for i, (date_str, entries) in enumerate(sorted_entries):
            # Add blank line before date (except for first entry)
            if i > 0:
                lines.append('')
            
            lines.append(f'- {date_str}:')
            
            # Group entries by category for better organization
            by_category = {}
            for entry in entries:
                cat = entry['category']
                if cat not in by_category:
                    by_category[cat] = []
                by_category[cat].append(entry)
            
            # Output in preferred category order
            category_order = ['feature', 'bugfix', 'change', 'maintenance', 'doc']
            
            for category in category_order:
                if category in by_category:
                    for entry in by_category[category]:
                        lines.append(f"  - [{entry['category']}] {entry['description']}")
                        for detail in entry['details']:
                            lines.append(f"    - {detail}")
        
        return '\n'.join(lines)
    
    def save_json_debug(self, data: Dict, filename: str = 'changelog_debug.json'):
        """Save parsed data to JSON for debugging"""
        with open(filename, 'w') as f:
            # Convert datetime objects to strings for JSON serialization
            json_data = data.copy()
            json.dump(json_data, f, indent=2, default=str)
        print(f"Debug data saved to {filename}")
    
    def restructure(self):
        """Main method to restructure the changelog"""
        print("Step 1: Parsing CHANGELOG.md to JSON...")
        data = self.parse_to_json()
        print(f"  Found {len(data['entries'])} unique dates")
        print(f"  Total raw entries: {len(data['raw_entries'])}")
        
        # Save debug JSON
        self.save_json_debug(data)
        
        print("\nStep 2: Cleaning duplicates...")
        cleaned_data = self.clean_entries(data)
        print(f"  Cleaned entries: {len(cleaned_data['entries'])} dates remain")
        
        print("\nStep 3: Sorting entries by date...")
        sorted_entries = self.sort_entries(cleaned_data)
        
        if sorted_entries:
            newest = sorted_entries[0][0]
            oldest = sorted_entries[-1][0]
            print(f"  Date range: {newest} → {oldest}")
        
        print("\nStep 4: Formatting back to markdown...")
        markdown = self.format_to_markdown(cleaned_data, sorted_entries)
        
        # Backup original
        backup_path = self.changelog_path.with_suffix('.md.backup')
        self.changelog_path.rename(backup_path)
        print(f"\nOriginal backed up to: {backup_path}")
        
        # Write new changelog
        self.changelog_path.write_text(markdown)
        print(f"✓ Restructured CHANGELOG saved to: {self.changelog_path}")
        
        # Summary statistics
        total_entries = sum(len(entries) for _, entries in sorted_entries)
        print(f"\nSummary:")
        print(f"  - Total dates: {len(sorted_entries)}")
        print(f"  - Total entries: {total_entries}")
        print(f"  - Average entries per date: {total_entries/len(sorted_entries):.1f}")


def main():
    restructurer = ChangelogRestructurer()
    restructurer.restructure()


if __name__ == '__main__':
    main()