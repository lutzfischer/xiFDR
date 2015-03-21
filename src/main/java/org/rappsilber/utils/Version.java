/*
 * Copyright 2015 Lutz Fischer <lfischer at staffmail.ed.ac.uk>.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package org.rappsilber.utils;

/**
 *
 * @author lfischer
 */
public class Version {
    /** the major number of the version*/
    public int major;
    /** the minor number of the version*/
    public int minor;
    /** build - svn-revision*/
    public int build;

    public Version(int major, int minor, int build) {
        this.major = major;
        this.minor = minor;
        this.build = build;
    }

    public Version(int major, int minor, String svn_refbuild) {
        this.major = major;
        this.minor = minor;
        this.build = Integer.parseInt(svn_refbuild.replaceAll("\\$Rev:\\s*", "").replaceAll("\\s*\\$", ""));
    }

    
    public String toString() {
        
         return major + "." + minor + "." + build;
    }

    public String toLongString() {
        
         return String.format("%d.%02d.%07d", major ,minor ,build);
    }
    
}
